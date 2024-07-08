#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# OVERVIEW
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Purpose: to compute the jSFS entries for simulated VCF datasets
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Run as: 
# python compute_jSFS.py --vcf [input vcf] --popfile [population manifest] --out [output filename prefix]
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------

#---------------------------
# IMPORT PACKAGES
#---------------------------
import allel
import warnings
warnings.simplefilter(action='ignore')
import pandas as pd
import numpy as np
import argparse
import subprocess
import gzip
import itertools
import os

#----------------------------------
# DEFINE THE MAIN SCRIPT FUNCTION
#----------------------------------

def main():

	#------------------------
	# Parse input arguments
	#------------------------

	print('Parsing command line arguments...')

	# Parse command-line arguments
	args = parse_args()
	# VCF to analyze
	vcf = args.vcf
	# Output file name to use to store results
	out = args.out
	# Popfile to use
	popfile = args.popfile

	print('Done.')

	print('Reading in the VCF file with scikit-allel...')

	# Read in the VCF file
	vcf_allel = allel.read_vcf(vcf, fields='*')
	# Get a list of samples in the VCF
	sample_vec = vcf_allel['samples']
	# Extract the genotype field
	# [NOTE] that genos[i] gives the genotypes of the ith variant in all individuals
	# [NOTE] that genos[i][j] gives the genotypes of the ith variant in individual j
	geno_vec = vcf_allel['calldata/GT']
	# Get the vector of chromosomes (1st col in VCF)
	chrom_vec = vcf_allel['variants/CHROM']
	# Get the vector of positions
	pos_vec = vcf_allel['variants/POS']
	# Get a list of the chromosomes represented in the VCF file
	chrom = np.unique(chrom_vec)

	# Print some output messages to ensure all is working properly
	print("Found samples:")
	print(sample_vec)
	print("Found " + str(len(pos_vec)) + " sites")
	print("Found chromosome:")
	print(chrom)

	print('Done')

	print('Creating population dictonary from popfile...')

	# Read in the popfile as a pandas CSV
	popfile = pd.read_csv(popfile, sep="\t", header=None)

	# Identify all unique populations represented in popfile, store names as list
	pop_ls = np.unique(popfile.iloc[:, 1].tolist())

	# Create a nested dictionary that will hold the population name as outermost key, then the sample names as nested inner keys
	pop_dict = {}

	# Loop through all populations represented in popfile to create dictionary
	for name in pop_ls:
		# Add an outer key for each population
		pop_dict[name] = {}
		# Within each population key, create a key for sample names, then populate accordingly
		pop_dict[name]['samples'] = popfile[popfile.iloc[:, 1]==name].iloc[:, 0].tolist()

		# Augment the popfile with additional info if running in summary statistic mode
		if(args.summary_stats=="yes"):
			# Also create a key for sample indices, then populate accordingly using their positions in the sample_vec
			pop_dict[name]['indices'] = [sample_vec.tolist().index(indv) for indv in pop_dict[name]['samples']]
			# Also create a key to store the corresponding subset of the geno_vec
			pop_dict[name]['geno_vec'] = geno_vec[:, pop_dict[name]['indices']]

	#----------------------------------------------------
	# Computing summary statistics for simulated outputs
	#----------------------------------------------------

	# Perform the jSFS calculation
	jSFS_df = jSFS(pop_ls, chrom, pop_dict, chrom_vec, pos_vec)

	# Add a column to the df that specifies the file name (containing relevant parameters)
	jSFS_df['params'] = out

	# Add a column to the df that specifies the chromosome name
	jSFS_df['chrom'] = chrom[0]

	# Write the full df to a tab-delimited file
	jSFS_df.to_string((out + '.jSFS'), index=False)

#---------------------
# DEFINE FUNCTIONS
#---------------------

def parse_args():
	"""
	Parse the command-line arguments
	"""
	# Call to argparse
	parser = argparse.ArgumentParser()

	# Passes VCF file to use (should be gzipped)
	parser.add_argument('--vcf', required=True)
	# Specifies name of output file
	parser.add_argument('--out', required=True)
	# Specifies popfile to use
	parser.add_argument('--popfile', required=True)

	# Parse and return arguments
	args = parser.parse_args()
	return args

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

def jSFS(pop_ls, chrom, pop_dict, chrom_vec, pos_vec):
	"""
	Compute the joint SFS for each pair of populations
	"""
	# Create a list to store the jSFS results for each population comparison
	jSFS_df_ls = []

	# Iterate through all pairwise comparisons of the populations specified in the popfile
	for i in range(0, len(pop_ls)):
		for j in range(i+1, len(pop_ls)):

			# Assign each population name to a variable
			pop1 = pop_ls[i]
			pop2 = pop_ls[j]

			# Create an empty df to store the jSFS results for this comparison
			jSFS_df = pd.DataFrame()

			print('Computing the jSFS for populations ' + pop1 + ' and ' + pop2 + '...')

			# Make the population 1 GenotypeArray object and the positions list for the given chromosome
			chrom_ga_pop1,chrom_pv_pop1 = make_ga(pop_dict[pop1]["geno_vec"], chrom_vec, pos_vec, chrom)
			# Make the population 2 GenotypeArray object and the positions list for the given chromosome
			chrom_ga_pop2,chrom_pv_pop2 = make_ga(pop_dict[pop2]["geno_vec"], chrom_vec, pos_vec, chrom)

			# Make an allele counts array for population 1 from the corresponding GenotypeArray object
			chrom_ac_pop1 = chrom_ga_pop1.count_alleles()
			# Make an allele counts array for population 2 from the corresponding GenotypeArray object
			chrom_ac_pop2 = chrom_ga_pop2.count_alleles()

			# Compute the jSFS as a 2D matric
			jSFS_matrix = allel.joint_sfs_folded(chrom_ac_pop1, chrom_ac_pop2)

			# Create a list to store pop1 minor allele counts
			pop1_mac = []
			# Create a list to store pop2 minor allele counts
			pop2_mac = []
			# Create a list to store corresponding number of SNPs in each bin
			num_snps = []

			# Loop through the matrix to format the output table
			for i in range(len(jSFS_matrix)):
				for j in range(len(jSFS_matrix[0])):
					# Append values to corresponging lists
					pop1_mac.append(i)
					pop2_mac.append(j)
					num_snps.append(jSFS_matrix[i][j])

			# Add these to the df
			jSFS_df['pop1_mac'] = pop1_mac
			jSFS_df['pop2_mac'] = pop2_mac
			jSFS_df['num_snps'] = num_snps

			# Create an ID for this comparison
			comp = pop1 + "-" + pop2

			# Add this info to the df
			jSFS_df['population'] = comp

			# Add this sub df to the list
			jSFS_df_ls.append(jSFS_df)

	# Concatenate all of the dfs together
	jSFS_df = pd.concat(jSFS_df_ls)

	return jSFS_df

#-------------------------
# RUN THE MAIN FUNCTION
#-------------------------

if __name__ == '__main__':
	main()








