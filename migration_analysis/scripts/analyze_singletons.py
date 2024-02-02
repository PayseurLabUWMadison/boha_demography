#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# OVERVIEW
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Purpose: to perform the singleton analysis on both simulated and empirical VCF datasets
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Run as: 
# python analyze_singletons.py --vcf [input vcf] --popfile [population manifest] --out [output filename prefix]
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------

#---------------------------
# IMPORT PACKAGES
#---------------------------
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

	print('Checking validity of sample and chromosome names...')

	# Need to replace all 0's in the sample names with something else (X here) so that PLINK won't throw an error

	# Variable to store output lines
	lines = []
	# Open the popfile
	with gzip.open(vcf, 'rt') as vcf_f:
		# Store lines into memory
		lines = vcf_f.readlines()
		# Iterate through the stored lines
		for num, line in enumerate(lines, 0):
			# Check to see if we're on the right header line
			if "#CHROM" in line:
				# Replace any 0's with X's
				lines[num] = line.replace("0", "X")
				# Remove any underscores
				lines[num] = lines[num].replace("_", "")
			# Check to see if we're onto the record lines
			elif ("#" not in line) and ("contig" not in line and "N" not in line) :
				# Append contig to the chrom ID (kinda sloppy)
				lines[num] = "contig" + line
	# Write the output popfile (overwrites the original)
	new_vcf = gzip.open(vcf, 'wt')
	new_vcf.writelines(lines)
	new_vcf.close()

	# Make a corresponding change to sample names in the popfile

	# Variable to store output lines
	lines = []
	# Open the popfile
	with open(popfile) as popfile_f:
		# Store lines into memory
		lines = popfile_f.readlines()
		# Iterate through the stored lines
		for num, line in enumerate(lines, 0):
			# Replace any 0's with X's
			lines[num] = line.replace("0", "X")
			# Remove any underscores
			lines[num] = lines[num].replace("_", "")
	# Write the output popfile (overwrites the original)
	new_popfile = open(popfile, 'w')
	new_popfile.writelines(lines)
	new_popfile.close()

	print('Done.')

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

	#--------------------------------------------------
	# Compute singleton statistics for each population
	#--------------------------------------------------
	
	# Create a list to store sub-dfs for the full singleton output
	singleton_df_ls = []

	# Loop through all populations to compute singletons
	for pop in pop_dict.keys():

		print('Identifying singletons in ' + pop + ' population...')

		# Create a "to-keep" list for this population that can be read by VCFtools
		keep = pd.DataFrame(list(zip(pop_dict[pop]['samples'], pop_dict[pop]['samples'])))
		keep_file = out + '.keep'
		keep.to_csv(keep_file, sep='\t', index=False, header=False)

		# Run the singleton function on the VCF dataset for the current population and capture the resulting dataframes
		s_df = singletons(vcf, keep_file, out)

		# Add a column to the sub-df that specifies population name
		s_df['population'] = pop
		# Add a column to the sub-df that specifies the file name (containing relevant parameters)
		s_df['params'] = out
		# Append the sub-df to the list
		singleton_df_ls.append(s_df)

	print('Done.')

	print('Writing output files...')

	# Concatenate all of the dfs contained in the list
	singleton_df = pd.concat(singleton_df_ls)

	# Write the full df to a tab-delimited file
	singleton_df.to_string((out + '.singletons'), index=False)

	print('Done.')

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
	parser.add_argument('--out', required=False)
	# Specifies popfile to use
	parser.add_argument('--popfile', required=False)

	# Parse and return arguments
	args = parser.parse_args()
	return args

def singletons(vcf, keep, out):
	"""
	Identify singletons in the VCF input

	returns: a dataframe
	"""
	print("Computing full singleton summary...")

	# Invoke the vcftools command
	subprocess.run(['sh', 'vcftools_singletons.sh', vcf, keep, out])
	# Capture the output of the vcftools ld command
	singleton_file = out + '.singletons'

	# Variable to store the output df
	s_df = ''

	# Reformat the singleton output

	# Check to see that the output file exists
	if(os.path.isfile(singleton_file)):
		# Read in the vcftools output, keeping only relevant columns
		s_df = pd.read_table(singleton_file, delim_whitespace=True, usecols=['CHROM', 'POS', 'SINGLETON/DOUBLETON', 'INDV'])
		# Rename the columns for clarity
		s_df.columns = ['chrom', 'pos', 'designation', 'indv']
	else:
		print('Issue with VCFtools command; not enough variants')

	return s_df

#-------------------------
# RUN THE MAIN FUNCTION
#-------------------------

if __name__ == '__main__':
	main()








