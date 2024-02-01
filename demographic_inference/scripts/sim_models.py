#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# OVERVIEW
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Purpose: Simulate sequence data under the best-fit demographic model
# Output: A formatted VCF file
# Sampling: Number of diploid individuals can be controlled with the nsamp options to match observed sample sizes
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Run as: 
# python sim_models.py --model [demographic model] --pop1_id [name of pop1] --pop1_nsamp [number of diploid individuals in pop1] --pop2_id [name of pop2] --pop2_nsamp [number of diploid individuals in pop2] --nsites [length of genomic element] --nreps [number of independent genomic elements] --recomb [per-bp, per-gen recombination rate] --mut [per-bp, per-gen recombination rate] 
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------

#---------------------------
# IMPORT PACKAGES
#---------------------------
import msprime
import tskit
import pandas as pd
from random import sample
import random
import numpy as np
import matplotlib.pyplot as pyplot
import argparse
import pylab

#----------------------------------
# DEFINE THE MAIN SCRIPT FUNCTION
#----------------------------------

def main():

	#-----------------
	# Set variables 
	#-----------------

	print("Reading in the data and intializing variables...")

	# Parse command-line arguments
	args = parse_args()

	# Set the demographic model for which we are simulating data
	model = args.model
	# Specify population names
	pop1_id = args.pop1_id
	pop2_id = args.pop2_id
	# Specify diploid sample sizes for each population
	pop1_nsamp = int(args.pop1_nsamp)
	pop2_nsamp = int(args.pop2_nsamp)
	# Set the length of the independent genomic segments to simulate
	nsites = int(args.nsites)
	# Set the number of independent genomic segments to simulate
	nreps = int(args.nreps)
	# Set the recombination rate
	recomb_rate = float(args.recomb)
	# Set the mutation rate
	mut_rate = float(args.mut)
	# These are the names of the output files
	vcf = model + "_mut" + str(args.mut) + "_rec" + str(args.recomb) + ".vcf"

	print("Done.")

	#---------------------------------
	# Construct the demographic model
	#---------------------------------

	print("Constructing demographic model...")

	# Construct the demographic model requested by user

	# Bumpkin-Peddocks 2D candidate
	elif (model == "simple_split_BP"):
		demography = simple_split_BP(pop1_id, pop2_id)
	# Bumpkin-WorldsEnd 2D candidate
	elif (model == "simple_split_BW"):
		demography = simple_split_BW(pop1_id, pop2_id)
	# Peddocks-WorldsEnd 2D candidate
	elif (model == "simple_split_PW"):
		demography = simple_split_PW(pop1_id, pop2_id)
	else:
		raise NameError('Model not defined!')

	print("Done.")

	#------------------------------------------------------------------------
	# Simulate the tree sequences, add mutations, and convert to a dataframe
	#------------------------------------------------------------------------

	print("Simulating sequence data...")
	# Use msprime to generate simulated data under above demographic model; output as VCF-like pandas dataframe
	df = sim_data(demography, pop1_id, pop2_id, pop1_nsamp, pop2_nsamp, nsites, nreps, recomb_rate, mut_rate)
	print("Done.")

	#-------------------------------
	# Writing VCF-like output
	#-------------------------------

	# Need a properly formatted header line
	print("Creating metadata header...")
	header = get_header(nreps, nsites)
	print("Done.")

	# Write full output to tab-delimited CSV file
	print("Writing VCF output...")
	with open(vcf, 'w') as outfile:
		outfile.write(header + "\n")
		df.to_csv(outfile, sep="\t", index=False)
	print("Done.")

#---------------------
# DEFINE FUNCTIONS
#---------------------

def parse_args():
	"""
	Parse the command-line arguments
	"""
	# Call to argparse
	parser = argparse.ArgumentParser()

	# Command-line variables

	# Determines which demographic model should be used
	parser.add_argument('--model', required=True)
	# Sets population names
	parser.add_argument('--pop1_id', required=True)
	parser.add_argument('--pop2_id', dest='pop2_id', default=False)
	# Determines how many diploid samples to simulate for each population
	parser.add_argument('--pop1_nsamp', required=True)
	parser.add_argument('--pop2_nsamp', dest='pop2_nsamp', default=False)	
	# Determines the length of each independent genomic segment
	parser.add_argument('--nsites', required=True)
	# Determines the number of genomic segments to simulate
	parser.add_argument('--nreps', required=True)
	# Determines the per-bp per-gen recombination rate
	parser.add_argument('--recomb', required=True)
	# Determines the per-bp per-gen mutation rate
	parser.add_argument('--mut', required=True)

	# Parse and return arguments
	args = parser.parse_args()
	return args

#-----------------------------------------
# Define demographic models
#-----------------------------------------

def simple_split_BP(pop1_id, pop2_id):
	"""
	simple_split_BP model
	return: an msprime demography object
	"""
	demography = msprime.Demography()
	demography.add_population(name=pop1_id, initial_size=7501)
	demography.add_population(name=pop2_id, initial_size=9465)
	demography.add_population(name="Ancestral", initial_size=470513)
	demography.add_population_split(time=7165, derived=[pop1_id, pop2_id], ancestral="Ancestral")

	return demography

def simple_split_BW(pop1_id, pop2_id):
	"""
	simple_split_BW model
	return: an msprime demography object
	"""
	demography = msprime.Demography()
	demography.add_population(name=pop1_id, initial_size=7257)
	demography.add_population(name=pop2_id, initial_size=67573)
	demography.add_population(name="Ancestral", initial_size=413589)
	demography.add_population_split(time=6911, derived=[pop1_id, pop2_id], ancestral="Ancestral")

	return demography

def simple_split_PW(pop1_id, pop2_id):
	"""
	simple_split_PW model
	return: an msprime demography object
	"""
	demography = msprime.Demography()
	demography.add_population(name=pop1_id, initial_size=9669)
	demography.add_population(name=pop2_id, initial_size=65483)
	demography.add_population(name="Ancestral", initial_size=464064)
	demography.add_population_split(time=7659, derived=[pop1_id, pop2_id], ancestral="Ancestral")

	return demography

#-----------------------------------------
# Simulation functions
#-----------------------------------------	

def sim_data(demography, pop1_id, pop2_id, pop1_nsamp, pop2_nsamp, nsites, nreps, recomb_rate, mut_rate):
	"""
	This function uses msprime to simulate the tree sequences, add mutations, and convert to a dataframe

	demography: this is the msprime demography object that was created
	nsamples: this is the number of diploid samples to simulates data for
	nsites: this is the length of each genomic segment to simulate
	nreps: this is the number of independent genomic segments to simulate
	recomb_rate: this is the per-bp per-gen recombination rate experienced by each genomic segment
	mut_rate: this is the per-bp per-gen mutation rate experienced by each genomic segment

	return: a pandas df storing the error-free simulated data in a VCF-like format
	"""
	# Create a list of lists (where each element is a row of the VCF output)
	data = list()

	# Loop through the required number of replicates, simulating the corresponding tree sequence objects and adding mutations to each
	for i in range(1, nreps+1):

		# Create the ts for the given replicate
		if(not pop2_id):
			ts = msprime.sim_ancestry(samples=pop1_nsamp, demography=demography, sequence_length=nsites, discrete_genome=True, recombination_rate=recomb_rate, ploidy=2)
		else:
			ts = msprime.sim_ancestry(samples={pop1_id: pop1_nsamp, pop2_id: pop2_nsamp}, demography=demography, sequence_length=nsites, discrete_genome=True, recombination_rate=recomb_rate, ploidy=2)
		# Then add mutations on this ts object
		mts = msprime.sim_mutations(ts, rate=mut_rate)
		# Store the current segment in VCF format
		vcf = mts.as_vcf(contig_id=i)
		# Turn the VCF output into a list of lines
		lines = vcf.split("\n")

		if(i % 100==0):
			print("Finished simulating replicate " + str(i))

		# Loop through the list of lines, breaking each into elements
		for line in lines:
			# Each element is a "cell" in the VCF
			elements = line.split("\t")

			# Only append rows that contain either the records or the column titles (we don't want the rest of the header info for now)
			# I also want to exclude rows/sites that have more than one derived allele (those will break the script)
			if (len(elements) > 1):

				# Count the number of ALT alleles (in case there are multi-allelic sites)
				if(len(elements[4].split(",")) < 2):

					# If this is the first replicate, output the header information along with VCF records to the data variable (so that we have proper column titles)
					if (i == 1):
						data.append(elements)
					
					# If this is not the first replicate, then make sure we don't include col titles since those have already been added
					elif (i > 1):
						if ("#CHROM" not in elements):
							data.append(elements)

	# Turn the data object that we generated for all replicates into a pandas df, using the first record as the column titles
	df = pd.DataFrame(data[1:], columns=data[0])

	return df

def get_header(nreps, nsites):
	"""
	This function generates an appropriate VCF-like header based on input parameters
	"""

	# Create a list to store the bookended and contig-specific header lines
	header_info = list()
	# Add the first several header elements to this list
	header_info.append("##fileformat=VCFv4.2")
	header_info.append("##source=tskit 0.5.4")
	header_info.append("##FILTER=<ID=PASS,Description=\"All filters passed\">")

	# For each independent chromosome, create a contig header line that specifies the length (given by nsites parameter)
	for i in range(0, nreps):
		# Construct VCF-like header string
		header_info.append("##contig=<ID=" + str(i+1) + ",length=" + str(nsites) + ">")

	# Add the last header element(s) to this list
	header_info.append("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">")

	# Join the header lines
	header = "\n".join(header_info)

	return header

#-------------------------
# RUN THE MAIN FUNCTION
#-------------------------

if __name__ == '__main__':
	main()




