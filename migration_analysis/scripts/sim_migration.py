#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# OVERVIEW
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Purpose: simulate data under best-fit two-population demographic models augmented with various migration scenarios
# Output: a single genomic segment with length specified by --nsites
# Sampling: want to recreate observed diploid sample sizes (Bumpkin-13, Peddocks-13, World's End-17) using --pop1_nsamp and --pop2_nsamp arguments
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Run as: 
# python sim_migration.py --model [demographic model] --mig_rate [symmetric migration rate with ghost population] --mig_time [proportion of divergence time over which migration occurs into present] --pop1_id [name of pop1] --pop2_id [name of pop2] --pop1_nsamp [diploid sample size of pop1] --pop2_nsamp [diploid sample size of pop2] --nsites [length of genomic element to simulate] --chrom [simulated chromosome ID] --recomb [per-bp, per-gen recombination rate] --mut [per-bp, per-gen mutation rate]
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------

#---------------------------
# IMPORT PACKAGES
#---------------------------
import msprime
import pandas as pd
import numpy as np
import argparse

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
	# Set the migration rate
	mig_rate = float(args.mig_rate)
	# Set the migration timing (as a fraction of total split time)
	mig_time = float(args.mig_time)

	# Specify population names
	pop1_id = args.pop1_id
	pop2_id = args.pop2_id
	# Specify diploid sample sizes for each population
	pop1_nsamp = int(args.pop1_nsamp)
	pop2_nsamp = int(args.pop2_nsamp)

	# Set the length of the independent genomic segments to simulate
	nsites = int(args.nsites)
	# This specifies the chromosome number that should be simulated
	chrom = int(args.chrom)
	# Set the recombination rate
	recomb_rate = float(args.recomb)
	# Set the mutation rate
	mut_rate = float(args.mut)

	# Output file name
	vcf = "Model-" + str(args.model) + "_MigRate-" + str(args.mig_rate) + "_MigTime-" + str(args.mig_time) + "_MutRate-" + str(args.mut) + "_RecRate-" + str(args.recomb) + "_Chrom-" + str(args.chrom) + ".vcf"

	print("Done.")

	#---------------------------------
	# Construct the demographic model
	#---------------------------------

	print("Constructing demographic model...")

	# Construct the demographic model requested by user
	if (model == "BP"):
		demography = BP(pop1_id, pop2_id, mig_rate, mig_time)
	elif (model == "BW"):
		demography = BW(pop1_id, pop2_id, mig_rate, mig_time)
	elif (model == "PW"):
		demography = PW(pop1_id, pop2_id, mig_rate, mig_time)
	else:
		raise NameError('Model not defined!')

	print("Done.")

	#------------------------------------------------------------------------
	# Simulate the tree sequences, add mutations, and convert to a dataframe
	#------------------------------------------------------------------------

	print("Simulating sequence data...")
	# Use msprime to generate simulated data under above demographic model; output as VCF-like pandas dataframe
	df = sim_data(demography, pop1_id, pop2_id, pop1_nsamp, pop2_nsamp, nsites, recomb_rate, mut_rate, chrom)
	print("Done.")

	#-------------------------------
	# Writing VCF-like output
	#-------------------------------

	# Need a properly formatted header line
	print("Creating metadata header...")
	header = get_header(chrom, nsites)
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
	# Determines what the symmetric migration rate parameter should be
	parser.add_argument('--mig_rate', required=True)
	# Determines when (as a fraction of total split time) migration should begin
	parser.add_argument('--mig_time', required=True)
	# Sets population names
	parser.add_argument('--pop1_id', required=True)
	parser.add_argument('--pop2_id', required=True)
	# Determines how many diploid samples to simulate for each population
	parser.add_argument('--pop1_nsamp', required=True)
	parser.add_argument('--pop2_nsamp', required=True)	
	# Determines the length of each independent genomic segment
	parser.add_argument('--nsites', required=True)
	# Determines which chromosome is being simulated
	parser.add_argument('--chrom', required=True)
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

def BP(pop1_id, pop2_id, mig_rate, mig_time):
	"""
	simple_split_BP model
	return: an msprime demography object
	"""
	# Define best-fit demographic parameters inferred with moments
	pop1_size = 7501
	pop2_size = 9465
	anc_size = 470513
	split_time = 7165

	# Initialize demography object
	demography = msprime.Demography()

	# Specify current sizes of each sampled population
	demography.add_population(name=pop1_id, initial_size=pop1_size)
	demography.add_population(name=pop2_id, initial_size=pop2_size)

	# Specify size of ancestral population
	demography.add_population(name="Ancestral", initial_size=anc_size)

	# Create a "ghost" population that has the same size as the ancestral population
	demography.add_population(name="Ghost", initial_size=anc_size)

	# Model population split, allowing the unsampled ghost population to diverge at same time as sampled population pair
	demography.add_population_split(time=split_time, derived=[pop1_id, pop2_id, "Ghost"], ancestral="Ancestral")

	# If migration parameters are specified, add migration to model
	if (mig_rate != 0) and (mig_time != 0):
		# Return message
		print("Modeling migration with rate = " + str(mig_rate) + " and time = " + str(mig_time) + " of total split time...")
		# Model the start and stop of symmetric migration between ghost population and pop1
		demography.add_symmetric_migration_rate_change(time=0, populations=[pop1_id,"Ghost"], rate=mig_rate)
		demography.add_symmetric_migration_rate_change(time=split_time*mig_time, populations=[pop1_id,"Ghost"], rate=0)
		# Model the start and stop of symmetric migration between ghost population and pop2
		demography.add_symmetric_migration_rate_change(time=0, populations=[pop2_id,"Ghost"], rate=mig_rate)
		demography.add_symmetric_migration_rate_change(time=split_time*mig_time, populations=[pop2_id,"Ghost"], rate=0)
	# If no migration parameters are specified, leave model as-is
	else:
		# Return message
		print("Modeling population divergence without migration...")

	# Output demography debug
	demography.sort_events()
	debug = demography.debug()
	print(debug)

	return demography

def BW(pop1_id, pop2_id, mig_rate, mig_time):
	"""
	simple_split_BW model
	return: an msprime demography object
	"""
	# Define best-fit demographic parameters inferred with moments
	pop1_size = 7257
	pop2_size = 67573
	anc_size = 413589
	split_time = 6911

	# Initialize demography object
	demography = msprime.Demography()

	# Specify current sizes of each sampled population
	demography.add_population(name=pop1_id, initial_size=pop1_size)
	demography.add_population(name=pop2_id, initial_size=pop2_size)

	# Specify size of ancestral population
	demography.add_population(name="Ancestral", initial_size=anc_size)

	# Create a "ghost" population that has the same size as the ancestral population
	demography.add_population(name="Ghost", initial_size=anc_size)

	# Model population split, allowing the unsampled ghost population to diverge at same time as sampled population pair
	demography.add_population_split(time=split_time, derived=[pop1_id, pop2_id, "Ghost"], ancestral="Ancestral")

	# If migration parameters are specified, add migration to model
	if (mig_rate != 0) and (mig_time != 0):
		# Return message
		print("Modeling migration with rate = " + str(mig_rate) + " and time = " + str(mig_time) + " of total split time...")
		# Model the start and stop of symmetric migration between ghost population and pop1
		demography.add_symmetric_migration_rate_change(time=0, populations=[pop1_id,"Ghost"], rate=mig_rate)
		demography.add_symmetric_migration_rate_change(time=split_time*mig_time, populations=[pop1_id,"Ghost"], rate=0)
		# Model the start and stop of symmetric migration between ghost population and pop2
		demography.add_symmetric_migration_rate_change(time=0, populations=[pop2_id,"Ghost"], rate=mig_rate)
		demography.add_symmetric_migration_rate_change(time=split_time*mig_time, populations=[pop2_id,"Ghost"], rate=0)
	# If no migration parameters are specified, leave model as-is
	else:
		# Return message
		print("Modeling population divergence without migration...")

	# Output demography debug
	demography.sort_events()
	debug = demography.debug()
	print(debug)

	return demography

def PW(pop1_id, pop2_id, mig_rate, mig_time):
	"""
	simple_split_PW model
	return: an msprime demography object
	"""
	# Define best-fit demographic parameters inferred with moments
	pop1_size = 9669
	pop2_size = 65483
	anc_size = 464064
	split_time = 7659

	# Initialize demography object
	demography = msprime.Demography()

	# Specify current sizes of each sampled population
	demography.add_population(name=pop1_id, initial_size=pop1_size)
	demography.add_population(name=pop2_id, initial_size=pop2_size)

	# Specify size of ancestral population
	demography.add_population(name="Ancestral", initial_size=anc_size)

	# Create a "ghost" population that has the same size as the ancestral population
	demography.add_population(name="Ghost", initial_size=anc_size)

	# Model population split, allowing the unsampled ghost population to diverge at same time as sampled population pair
	demography.add_population_split(time=split_time, derived=[pop1_id, pop2_id, "Ghost"], ancestral="Ancestral")

	# If migration parameters are specified, add migration to model
	if (mig_rate != 0) and (mig_time != 0):
		# Return message
		print("Modeling migration with rate = " + str(mig_rate) + " and time = " + str(mig_time) + " of total split time...")
		# Model the start and stop of symmetric migration between ghost population and pop1
		demography.add_symmetric_migration_rate_change(time=0, populations=[pop1_id,"Ghost"], rate=mig_rate)
		demography.add_symmetric_migration_rate_change(time=split_time*mig_time, populations=[pop1_id,"Ghost"], rate=0)
		# Model the start and stop of symmetric migration between ghost population and pop2
		demography.add_symmetric_migration_rate_change(time=0, populations=[pop2_id,"Ghost"], rate=mig_rate)
		demography.add_symmetric_migration_rate_change(time=split_time*mig_time, populations=[pop2_id,"Ghost"], rate=0)
	# If no migration parameters are specified, leave model as-is
	else:
		# Return message
		print("Modeling population divergence without migration...")

	# Output demography debug
	demography.sort_events()
	debug = demography.debug()
	print(debug)

	return demography

#-----------------------------------------
# Simulation functions
#-----------------------------------------	

def sim_data(demography, pop1_id, pop2_id, pop1_nsamp, pop2_nsamp, nsites, recomb_rate, mut_rate, chrom):
	"""
	This function uses msprime to simulate the tree sequences, add mutations, and convert to a dataframe

	demography: this is the msprime demography object that was created
	nsamples: this is the number of diploid samples to simulates data for
	nsites: this is the length of each genomic segment to simulate
	recomb_rate: this is the per-bp per-gen recombination rate experienced by each genomic segment
	mut_rate: this is the per-bp per-gen mutation rate experienced by each genomic segment

	return: a pandas df storing the simulated data in a VCF-like format
	"""
	# Create a list of lists (where each element is a row of the VCF output)
	data = list()

	# Simulate the tree sequence object and add mutations

	# Create the ts for the given replicate
	ts = msprime.sim_ancestry(samples={pop1_id: pop1_nsamp, pop2_id: pop2_nsamp}, demography=demography, sequence_length=nsites, discrete_genome=True, recombination_rate=recomb_rate, ploidy=2)
	# Then add mutations on this ts object
	mts = msprime.sim_mutations(ts, rate=mut_rate)
	# Store the current segment in VCF format
	vcf = mts.as_vcf(contig_id=chrom)
	# Turn the VCF output into a list of lines
	lines = vcf.split("\n")

	print("Finished simulating chromosome " + str(chrom))

	# Loop through the list of lines, breaking each into elements
	for line in lines:
		# Each element is a "cell" in the VCF
		elements = line.split("\t")

		# Only append rows that contain either the records or the column titles (we don't want the rest of the header info for now)
		# I also want to exclude rows/sites that have more than one derived allele (those will break the script)
		if (len(elements) > 1):

			# Count the number of ALT alleles (in case there are multi-allelic sites)
			if(len(elements[4].split(",")) < 2):

				# Write all lines (header plus VCF records)
				data.append(elements)

	# Turn the data object that we generated for all replicates into a pandas df, using the first record as the column titles
	df = pd.DataFrame(data[1:], columns=data[0])

	return df

def get_header(chrom, nsites):
	"""
	This function generates an appropriate VCF-like header based on input parameters
	"""

	# Create a list to store the bookended and contig-specific header lines
	header_info = list()
	# Add the first several header elements to this list
	header_info.append("##fileformat=VCFv4.2")
	header_info.append("##source=tskit 0.5.4")
	header_info.append("##FILTER=<ID=PASS,Description=\"All filters passed\">")

	# For simulated chromosome, create a contig header line that specifies the length (given by nsites parameter)
	header_info.append("##contig=<ID=" + str(chrom) + ",length=" + str(nsites) + ">")

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




