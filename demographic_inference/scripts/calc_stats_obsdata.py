#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# OVERVIEW
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Purpose: Computes summary statistics in 5kbp non-overlapping windows from observed sequence data 
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Run as: 
# python calc_stats_obsdata.py --vcf [VCF for observed sequence data] --mask [corresponding mask file of inaccessible sites/regions] --out [output prefix] --[pi, tajD, or fst]
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------

#---------------------------
# IMPORT PACKAGES
#---------------------------
import allel
import pandas as pd
import numpy as np
import argparse
from scipy.spatial.distance import squareform
from scipy.spatial.distance import pdist
import subprocess

#----------------------------------
# DEFINE THE MAIN SCRIPT FUNCTION
#----------------------------------

def main():

	#-----------------
	# Set variables 
	#-----------------

	print("--------------------------------------------------------------")
	print("Analyzing inputs...")
	print("--------------------------------------------------------------")

	print("Intializing variables...")

	# Parse command-line arguments
	args = parse_args()
	# VCF to analyze
	invcf = args.vcf
	# Mask file to use
	inmask = args.mask
	# Output file prefix
	out = args.out

	print("Reading in the VCF file...")

	# Read in the VCF file
	vcf = allel.read_vcf(invcf, fields='*')
	# Get a list of samples in the VCF
	sample_vec = vcf['samples']
	# Extract the genotype field
	# [NOTE] that genos[i] gives the genotypes of the ith variant in all individuals
	# [NOTE] that genos[i][j] gives the genotypes of the ith variant in individual j
	geno_vec = vcf['calldata/GT']
	# Get the vector of chromosomes (1st col in VCF)
	chrom_vec = vcf['variants/CHROM']
	# Get the vector of positions
	pos_vec = vcf['variants/POS']
	# Get a list of the chromosomes represented in the VCF file
	chrom_ls = np.unique(chrom_vec)

	print("Found samples:")
	print(sample_vec)
	print("Found " + str(len(pos_vec)) + " sites")
	print("Found chromosomes:")
	print(chrom_ls)

	print("Parsing mask file...")

	# Parse the masked file to create a dictionary of start and stop indices for the lines pertaining to each contig
	contig_dict = parse_mask(inmask)

	#----------------------------------------------------------------------
	# Determine which summary statistics should be run based on user input
	#----------------------------------------------------------------------

	# Check to see whether popfile exists
	if(args.popfile is None):

		print("No popfile detected; running in single-population mode")

		# Perform each single population summary statistic calculation specified on cmd line
		if(args.pi):
			pi(out, chrom_ls, contig_dict, inmask, geno_vec, chrom_vec, pos_vec)
		if(args.tajD):
			tajD(out, chrom_ls, geno_vec, chrom_vec, pos_vec)

	# If popfile specified, run in multipop mode
	elif(args.popfile is not None):
		
		print("Popfile detected; running in multi-population mode")
		print("Creating population dictonary from popfile...")

		# Read in the popfile
		popfile = pd.read_csv(args.popfile, sep="\t", header=None)

		# Identify all unique populations represented in file, store names as list
		pop_ls = np.unique(popfile.iloc[:, 1].tolist())

		# Create a nested dictionary that will hold the population name as outermost key, then the sample names and corresponding indices in sample_vec as nested inner keys
		pop_dict = {}

		# Loop through all populations represented in popfile to create dictionary
		for name in pop_ls:
			# Add an outer key for each population
			pop_dict[name] = {}
			# Within each population key, create a key for sample names, then populate accordingly
			pop_dict[name]["samples"] = popfile[popfile.iloc[:, 1]==name].iloc[:, 0].tolist()
			# Also create a key for sample indices, then populate accordingly using their positions in the sample_vec
			pop_dict[name]["indices"] = [sample_vec.tolist().index(indv) for indv in pop_dict[name]["samples"]]
			# Also create a key to store the corresponding subset of the geno_vec
			pop_dict[name]["geno_vec"] = geno_vec[:, pop_dict[name]["indices"]]

		# Perform each multi-population summary statistic calculation specified on cmd line
		if(args.fst):
			fst(pop_ls, out, chrom_ls, pop_dict, chrom_vec, pos_vec)

#---------------------
# DEFINE FUNCTIONS
#---------------------

def parse_args():
	"""
	Parse the command-line arguments
	"""
	# Call to argparse
	parser = argparse.ArgumentParser()

	# File inputs

	# Passes VCF file to use (should be gzipped)
	parser.add_argument('--vcf', required=True)
	# Passes FASTA like mask file (mask character should be 'N', file should not be gzipped)
	parser.add_argument('--mask', required=True)
	# Passes the popfile to use (identical to dadi popfile format)
	parser.add_argument('--popfile', required=False)

	# Program options

	# Single population nucleotide diversity
	parser.add_argument('--pi', dest='pi', default=False, action='store_true')
	# Single population Tajima's D
	parser.add_argument('--tajD', dest='tajD', default=False, action='store_true')
	# Multi population fst 
	parser.add_argument('--fst', dest='fst', default=False, action='store_true')

	# Output files

	# Passes the prefix to use for output
	parser.add_argument('--out', required=True)

	# Parse and return arguments
	args = parser.parse_args()
	return args

def parse_mask(inmask):
	"""
	Parse the masked file to create a dictionary of start and stop indices for the lines pertaining to each contig
	"""
	# Stores a list of contigs found in the mask file
	contig = []
	# Stores a list of start line indices for each contig found in mask file ([NOTE] this is index of header line ">")
	start_idx = []
	# Stores a list of stop line indices for each contig found in mask file
	stop_idx = []

	# Get contig name, start, and stop positions to use in further parsing the mask file
	with open(inmask) as m: 
		# Loop through each line in mask file, track line content and index
		for idx, line in enumerate(m.readlines()):
			# Check that line is a header line
			if ">" in line:
				# Strip extraneous info to get contig name
				line = line.strip('\n')
				line = line.strip('>')
				# Append contig name to list
				contig.append(line)
				# Append the line number this header is found on to start_idx list
				start_idx.append(idx)
				# If this isn't the first header in the file, then append the previous position to the stop_idx list
				if idx > 0:
					stop_idx.append(idx - 1)
		# Once we've reached the end of the file, append the index of the last line to the stop_idx list
		stop_idx.append(idx)

	# Create a dictionary that contains each contig as the key and the [start, stop] positions as the value
	contig_dict = dict(zip(contig, [list(tup) for tup in zip(start_idx, stop_idx)]))
	return contig_dict

def make_ga(geno_vec, chrom_vec, pos_vec, this_chrom):
	"""
	Make the GenotypeArray object and the positions list for the given chromosome
	"""
	# This will store the genotype array for the given chromosome
	chrom_ga = []
	# This will store the positions vector for the given chromosome
	chrom_pv = []

	# Enumerate through genotype field, storing index (i) and genotypes (g) at the ith variant in all individuals
	for i, g in enumerate(geno_vec):
		# If this variant occurs on the chromosome of interest...
		if chrom_vec[i]==this_chrom:
			# Take the genotype of *all* individuals at this site
			chrom_ga.append(g)
			# Take the position for the variant
			chrom_pv.append(pos_vec[i])

	# Turn the collection of genotypes on the given chromosome into a GA object
	chrom_ga = allel.GenotypeArray(chrom_ga)
	return chrom_ga,chrom_pv

def make_accessible(contig_dict, this_chrom, inmask):
	"""
	For *only* the given chromosome, create an is_accessible array
	"""
	# Store the start/stop indices for the chromosome of interest using the key name of the chromosome
	coords = contig_dict[this_chrom]

	# Create a boolean array that stores whether the site is masked (False) or whether it is observed (True)
	accessible = []

	# Re-open the mask file
	with open(inmask) as m:
		# Only look at the lines pertaining to the given chromsome of interest (indices are a bit funky here- might want to double-check this)
		lines = m.readlines()[(coords[0] + 1):(coords[1] + 1)]
		# Loop through the lines in the mask file pertaining to the sites on the given chromosome
		for line in lines:
			# Turn each line of nucleotides to a list of elements
			line = line.strip('\n')
			elements = list(line)
			# Convert the list of elements to a boolean array denoting the mask status
			boolements = [(e != 'N') for e in elements]
			# Unlist the boolean array to append the observations back to the accessible array
			[accessible.append(b) for b in boolements]
	return accessible

def pi(out, chrom_ls, contig_dict, inmask, geno_vec, chrom_vec, pos_vec):
	""" 
	For all available chromsomes, compute windowed pi
	"""
	print("--------------------------------------------------------------")
	print("Beginning nucleotide diversity analysis...")
	print("--------------------------------------------------------------")

	# Output name for the pi file
	pi_file = out + ".pi"

	# Create a list to hold individual pandas df that contain the window-based pi estimates for each chrom
	pi_df_ls = []

	# Loop through all found chromosomes
	for chrom_idx in range(0, len(chrom_ls)):

		# Name of current chromosome
		this_chrom = chrom_ls[chrom_idx]

		print("Analyzing chromosome " + this_chrom + "...")

		# For *only* the given chromosome, create an is_accessible array
		accessible = make_accessible(contig_dict, this_chrom, inmask)

		# Make the GenotypeArray object and the positions list for the given chromosome
		chrom_ga,chrom_pv = make_ga(geno_vec, chrom_vec, pos_vec, this_chrom)

		# Make an allele counts array from the GenotypeArray object
		chrom_ac = chrom_ga.count_alleles()
		
		# Make sure there is at least one variant site in this chromosome
		if(chrom_ga.n_variants > 1):
			# Compute windowed nucleotide diversity, creating relevant df columns
			pi,windows,n_bases,counts = allel.windowed_diversity(chrom_pv, chrom_ac, size=5000, is_accessible=accessible, fill="NA")

			# Extract start and stop positions from windows output
			start = np.split(windows, 2, axis=1)[0].flatten()
			stop = np.split(windows, 2, axis=1)[1].flatten()

			# Create pandas df to hold the pi results
			pi_df = pd.DataFrame({'chrom': this_chrom, 'pi': pi, 'start': start, 'stop': stop, 'n_bases': n_bases, 'counts': counts})
			
			# Append the window-based pi estimates from this chromosome to the list
			pi_df_ls.append(pi_df)
		else:
			print("Not enough variants on chromosome " + this_chrom)

	print("Compiling results...")

	# Concatenate all of the dfs contained in the pi_df_ls list
	pi_df = pd.concat(pi_df_ls)

	print("Writing results to files...")

	# Write the pi df to a tab-delimited CSV
	pi_df.to_csv(pi_file, sep="\t", index=False)

def tajD(out, chrom_ls, geno_vec, chrom_vec, pos_vec):
	""" 
	For all available chromsomes, compute windowed Tajima's D
	"""
	print("--------------------------------------------------------------")
	print("Beginning Tajima's D analysis...")
	print("--------------------------------------------------------------")

	# Output name for the tajD file
	tajD_file = out + ".tajD"

	# Create a list to hold individual pandas df that contain the window-based tajD estimates for each chrom
	tajD_df_ls = []

	# Loop through all found chromosomes
	for chrom_idx in range(0, len(chrom_ls)):

		# Name of current chromosome
		this_chrom = chrom_ls[chrom_idx]

		print("Analyzing chromosome " + this_chrom + "...")

		# Make the GenotypeArray object and the positions list for the given chromosome
		chrom_ga,chrom_pv = make_ga(geno_vec, chrom_vec, pos_vec, this_chrom)

		# Make an allele counts array from the GenotypeArray object
		chrom_ac = chrom_ga.count_alleles()

		# Make sure there is at least one variant site in this chromosome
		if(chrom_ga.n_variants > 1):
			# Compute windowed Tajima's D, creating relevant df columns
			tajD,windows,counts = allel.windowed_tajima_d(chrom_pv, chrom_ac, size=5000)

			# Extract start and stop positions from windows output
			start = np.split(windows, 2, axis=1)[0].flatten()
			stop = np.split(windows, 2, axis=1)[1].flatten()

			# Create pandas df to hold the tajD results
			tajD_df = pd.DataFrame({'chrom': this_chrom, 'tajD': tajD, 'start': start, 'stop': stop, 'counts': counts})
			
			# Append the window-based tajD estimates from this chromosome to the list
			tajD_df_ls.append(tajD_df)
		else:
			print("Not enough variants on chromosome " + this_chrom)

	print("Compiling results...")

	# Concatenate all of the dfs contained in the tajD_df_ls list
	tajD_df = pd.concat(tajD_df_ls)

	print("Writing results to files...")

	# Write the tajD df to a tab-delimited CSV
	tajD_df.to_csv(tajD_file, sep="\t", index=False)

def fst(pop_ls, out, chrom_ls, pop_dict, chrom_vec, pos_vec):
	"""
	Compute Hudson's Fst in windows for each pair of populations for all chromosomes
	"""
	print("--------------------------------------------------------------")
	print("Beginning Fst analysis...")
	print("--------------------------------------------------------------")

	# Iterate through all pairwise comparisons of the populations specified in the popfile
	for i in range(0, len(pop_ls)):
		for j in range(i+1, len(pop_ls)):

			# Assign each population name to a variable
			pop1 = pop_ls[i]
			pop2 = pop_ls[j]

			# Construct the filename for this comparison
			fst_file = pop1 + "_" + pop2 + "_" + out + ".fst"

			# Create a list to hold the dfs produced from this comparison for each chromosome
			fst_df_ls = []

			print("Comparing populations: " + pop1 + " and " + pop2)

			# Loop through all found chromosomes
			for chrom_idx in range(0, len(chrom_ls)):

				# Name of current chromosome
				this_chrom = chrom_ls[chrom_idx]

				print("Analyzing chromosome " + this_chrom + "...")

				# Make the population 1 GenotypeArray object and the positions list for the given chromosome
				chrom_ga_pop1,chrom_pv_pop1 = make_ga(pop_dict[pop1]["geno_vec"], chrom_vec, pos_vec, this_chrom)
				# Make the population 2 GenotypeArray object and the positions list for the given chromosome
				chrom_ga_pop2,chrom_pv_pop2 = make_ga(pop_dict[pop2]["geno_vec"], chrom_vec, pos_vec, this_chrom)

				# Make an allele counts array for population 1 from the corresponding GenotypeArray object
				chrom_ac_pop1 = chrom_ga_pop1.count_alleles()
				# Make an allele counts array for population 2 from the corresponding GenotypeArray object
				chrom_ac_pop2 = chrom_ga_pop2.count_alleles()

				# Compute Hudsons Fst, creating relevant df columns
				fst, windows, counts = allel.windowed_hudson_fst(chrom_pv_pop1, chrom_ac_pop1, chrom_ac_pop2, size=5000, fill="NA")

				# Extract start and stop positions from windows output
				start = np.split(windows, 2, axis=1)[0].flatten()
				stop = np.split(windows, 2, axis=1)[1].flatten()

				# Create pandas df to hold the Fst results
				fst_df = pd.DataFrame({'chrom': this_chrom, 'Fst': fst, 'start': start, 'stop': stop, 'counts': counts})
				
				# Append this df to the Fst df list
				fst_df_ls.append(fst_df)

			print("Compiling results...")

			# Concatenate all of the dfs contained in the fst_df_ls list
			fst_df = pd.concat(fst_df_ls)

			print("Writing results to files...")

			# Write the fst df to a tab-delimited CSV
			fst_df.to_csv(fst_file, sep="\t", index=False)

#-------------------------
# RUN THE MAIN FUNCTION
#-------------------------

if __name__ == '__main__':
	main()




