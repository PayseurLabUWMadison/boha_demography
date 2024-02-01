#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# OVERVIEW
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Purpose: estimate parameter uncertainties for specific cases (semi-hardcoded)
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Run as: 
# python est_uncert.py --vcf [filtered VCF used for original inference] --popfile [population manifest] --model [model] --popt [param1 param2 ... paramN]
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------

#---------------------------
# IMPORT PACKAGES
#---------------------------
import pandas as pd
from random import sample
import moments
import dill as pickle 
import random
import numpy as np
import matplotlib.pyplot as pyplot
import argparse
import pylab
import glob

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

	# Variable for holding the name of the input VCF file
	vcf = args.vcf
	# Variable for holding input popfile
	popfile = args.popfile
	# Create a variable to hold the model ID
	model_id = args.model
	# Create a variable to hold the parameter estimates
	popt_est = args.popt
	# Create a variable to hold the parameter estimates as a list of floats
	popt = [float(p) for p in popt_est]

	print("Done.")

	#---------------------------------------------
	# Constructing bootstrapped frequency spectra
	#---------------------------------------------

	# Run function to create bootstrapped FS from input VCF
	all_boot, original_fs = create_boots(vcf, popfile)

	#---------------------------------
	# Compute parameter uncertainties
	#---------------------------------

	# Determine which demographic model we're estimating uncertainties for
	if (model_id=="split"):
		model_func = split
	else:
		print("Error: no such model!")

	print("Original parameter estimates:")
	print(popt)

	print("Original theta:")
	# Determine the sample sizes from the original FS
	ns = original_fs.sample_sizes
	# Calculate the best-fit model FS
	model = model_func(popt, ns)
	# Calculate the corresponding theta value
	theta = moments.Inference.optimal_sfs_scaling(model, original_fs)
	print(theta)

	print("Estimating parameter uncertainties...")

	# Loop through different possible step sizes to ensure uncertainty estimates are stable
	for eps in [0.01, 0.001, 0.0001]:

		print("#------------------------------------------------------------")
		print("Using step size = " + str(eps))

		# Function returns standard deviation of parameter values (including theta– listed last) along with the full GIM to use in propogating uncertainties
		uncert, GIM = moments.Godambe.GIM_uncert(model_func, all_boot, popt, original_fs, log=False, multinom=True, eps=eps, return_GIM=True)

		print("Standard deviations of parameter values (plus theta):")
		print(uncert)

		print("Full GIM matrix:")
		print(GIM)

		print("Inverse of GIM matrix:")
		print(np.linalg.inv(GIM))
		print("#------------------------------------------------------------")


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

	# Determines which VCF to use
	parser.add_argument('--vcf', required=True)
	# Determines which popfile to use
	parser.add_argument('--popfile', required=True)
	# This is the demographic model we wish to fit 
	parser.add_argument("--model")
	# These are the best-fit parameter estimates (excluding theta and the ll) for the given model
	parser.add_argument("--popt", default=[], nargs="+")

	# Parse and return arguments
	args = parser.parse_args()
	return args

def create_boots(vcf, popfile):
	"""
	This function takes in a VCF dataset and creates bootstrap replicates to use in LRT
	"""
	# Read in the provided popfile as a pandas df
	popf_df = pd.read_csv(popfile, sep='\t', header=None)

	# Identify all unique populations in the popfile
	pop_ls = np.unique(popf_df.iloc[:,1].tolist())
	print("Populations detected:")
	print(pop_ls)

	# Construct a list for the samples that correspond to each unique population
	samp_ls = [popf_df[popf_df[1] == p][0].tolist() for p in pop_ls]
	print("Corresponding samples:")
	print(samp_ls)

	# Construct a list for the sample sizes to use in the projection argument
	proj_ls = [2*len(s) for s in samp_ls]
	print("Corresponding sample sizes:")
	print(proj_ls)

	print("Constructing data dictionary from input VCF...")

	# Read in the input vcf and construct moments dictionary object
	dd = moments.Misc.make_data_dict_vcf(vcf, popfile)

	print("Generating bootstraps...")

	# Create folded bootstrap replicats
	# And now write the output to a directory
	boot_dir = popfile + "_uncert_boots"
	moments.Misc.bootstrap(dd, pop_ls, proj_ls, polarized = False, num_boots=100, save_dir=boot_dir)

	# Load the saved bootstrap replicates
	boots_fids = glob.glob(boot_dir+'/*')
	all_boot = [moments.Spectrum.from_file(fid) for fid in boots_fids]

	print("Bootstraps generated.")

	print("Generating full frequency spectrum...")

	original_fs = moments.Spectrum.from_data_dict(dd, pop_ls, projections = proj_ls, polarized = False)

	print("Full frequency spectrum generated.")

	return all_boot, original_fs

# Define the demographic function for a simple population split
def split(params, ns):
	"""
	Parameter values:
	1. nu1 = the ratio of population 1 size to the ancestral population size
	2. nu2 = the ratio of population 2 size to the ancestral population size
	3. T = the time since the population split
	Overview: the ancestral population splits into population 1 with size nu1 and population 2 with size nu2 T generations from the present.

	[NOTE] all parameters are given relative to the reference size. The ancestral reference size must be calculated from theta
	"""
	# Initialize the parameter values
	nu1, nu2, T = params
	# Create an intital 1D AFS that represents our equilibrium ancestral population (use combined sample sizes)
	sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
	# Turn it into a FS object
	fs = moments.Spectrum(sts)
	# Create a 2D FS
	fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
	# After the splitting T generations ago, the populations take constant sizes nu1 and nu2
	fs.integrate([nu1, nu2], T)
	# Return the frequency spectrum obtained from the model
	return fs

#-------------------------
# RUN THE MAIN FUNCTION
#-------------------------

if __name__ == '__main__':
	main()




