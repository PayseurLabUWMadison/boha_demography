#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# OVERVIEW
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Purpose: perform LRT for specific cases (hardcoded)
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Run as: 
# python model_lrt.py --vcf [filtered VCF used for original inference] --popfile [population manifest]
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

	print("Done.")

	#---------------------------------------------
	# Constructing bootstrapped frequency spectra
	#---------------------------------------------

	# Run function to create bootstrapped FS from input VCF
	all_boot, original_fs = create_boots(vcf, popfile)

	#-----------------------------------
	# Set-up the likelihood ratio tests
	#-----------------------------------

	print("Setting up likelihood ratio tests...")

	# Determine the sample sizes from the original FS
	ns = original_fs.sample_sizes

	# Specify the best-fit parameter values for each model
	popt_split = np.array([0.0208344208194941, 0.141107010682326, 0.00825242250464207]) 
	popt_split_mig = np.array([0.021196261553823, 0.142776691439273, 0.00845189815028967, 0.303084939458155])
	popt_split_exp_const = np.array([1.34809211617666, 0.00581068160712739, 0.19729596102029, 0.0124144912526951])
	popt_split_exp_const_mig = np.array([2.51767553586668, 0.00667600922587468, 0.225385923487663, 0.0159919224923917, 2.52965012938893])

	# Calculate the best-fit model FS for each contrasting model
	model_fs_split = split(popt_split, ns)
	model_fs_split_mig = split_mig(popt_split_mig, ns)
	model_fs_split_exp_const = split_exp_const(popt_split_exp_const, ns)
	model_fs_split_exp_const_mig = split_exp_const_mig(popt_split_exp_const_mig, ns)

	# Calculate the likelihood of the data given each model FS
	ll_split = moments.Inference.ll_multinom(model_fs_split, original_fs)
	ll_split_mig = moments.Inference.ll_multinom(model_fs_split_mig, original_fs)
	ll_split_exp_const = moments.Inference.ll_multinom(model_fs_split_exp_const, original_fs)
	ll_split_exp_const_mig = moments.Inference.ll_multinom(model_fs_split_exp_const_mig, original_fs)

	print("Ready to perform LRT.")

	#------------------------------------
	# Perform the likelihood ratio tests
	#------------------------------------

	# Comparing simple split (with and without migration)
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	print("#-------------------------------------------")
	print("Comparing simple split (with and without migration)")

	# Use simplified model parameters as input into more complex model
	# In this case, we insert 0 for the migration rate parameter (which comes last in the popt list)
	p_lrt = np.append(popt_split, 0)

	print("Inputting simpler model estimates into more complex model:")
	print(p_lrt)

	# Determine which of the indices is nested (should be last one)
	#nested_index = p_lrt.size - 1
	nested_index = 3

	# Need to tell the LRT function which indices correspond to parameters absent in simplified model
	adj = moments.Godambe.LRT_adjust(split_mig, all_boot, p_lrt, original_fs, nested_indices=[nested_index], multinom=True)
	print('Adjustment value: ' + str(adj))

	# Multiply the adjustment by the difference in likelihoods between more complex model - less complex model
	D_adj = adj*2*(ll_split_mig - ll_split)
	print('Adjusted D statistic: ' + str(D_adj))

	# Calculate p-value from chi^2 distribution
	pval = moments.Godambe.sum_chi2_ppf(D_adj, weights=(0.5,0.5))
	print('p-value for rejecting simpler no-migration model: ' + str(pval))

	# Comparing split_exp_const (with and without migration)
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	print("#-------------------------------------------")
	print("Comparing split_exp_const (with and without migration)")

	# Use simplified model parameters as input into more complex model
	# In this case, we insert 0 for the migration rate parameter (which comes last in the popt list)
	p_lrt = np.append(popt_split_exp_const, 0)

	print("Inputting simpler model estimates into more complex model:")
	print(p_lrt)

	# Determine which of the indices is nested (should be last one)
	#nested_index = p_lrt.size - 1
	nested_index = 4

	# Need to tell the LRT function which indices correspond to parameters absent in simplified model
	adj = moments.Godambe.LRT_adjust(split_exp_const_mig, all_boot, p_lrt, original_fs, nested_indices=[nested_index], multinom=True)
	print('Adjustment value: ' + str(adj))

	# Multiply the adjustment by the difference in likelihoods between more complex model - less complex model
	D_adj = adj*2*(ll_split_mig - ll_split)
	print('Adjusted D statistic: ' + str(D_adj))

	# Calculate p-value from chi^2 distribution
	pval = moments.Godambe.sum_chi2_ppf(D_adj, weights=(0.5,0.5))
	print('p-value for rejecting simpler no-migration model: ' + str(pval))

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

	# Create folded bootstrap replicatse
	# And now write the output to a directory
	# moments.Misc.bootstrap(dd, pop_ls, proj_ls, polarized = False, num_boots=100, save_dir="peddocks_worlds_end_lrt_boots")

	# Load the saved bootstrap replicates
	boots_fids = glob.glob('peddocks_worlds_end_lrt_boots/*')
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

# Define the demographic function for a population split followed by exponential change in population 1 and no change in population 2
def split_exp_const(params, ns):
	"""
	Parameter values:
	1. nu1B = the ratio of population 1's bottleneck size to the ancestral population size
	2. nu1C = the ratio of population 1's current size to the ancestral population size
	3. nu2C = the ratio of population 2's current size to the ancestral population size
	4. T = the time since the population split
	Overview: the ancestral population splits into population 1 with size nu1B and population 2 with size nu2C T generations from the present.
		After the split, the population 1 undergoes exponential change to nu1C and population stays constant at nu2C.

	[NOTE] all parameters are given relative to the reference size. The ancestral reference size must be calculated from theta
	"""
	# Initialize the parameter values
	nu1B, nu1C, nu2C, T = params
	# Create an intital 1D AFS that represents our equilibrium ancestral population (use combined sample sizes)
	sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
	# Turn it into a FS object
	fs = moments.Spectrum(sts)
	# Create a 2D FS
	fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
	# Create a function that describes exponential growth in pop 1 from size nu1B to nu1C over T generations
	nu1_func = lambda t: nu1B * (nu1C/nu1B) ** (t/T)
	# Make a general nu_func
	nu_func = lambda t: [nu1_func(t), nu2C]
	# After the splitting T generations ago, pop 1 changes exponentially according to function while pop 2 stays the same
	fs.integrate(nu_func, T)
	# Return the frequency spectrum obtained from the model
	return fs

# Define the demographic function for a population split followed by constant sizes for both population 1 and 2 with symmetric migration
def split_mig(params, ns):
	"""
	Parameter values:
	1. nu1 = the ratio of population 1's size to the ancestral population size
	2. nu2 = the ratio of population 2's size to the ancestral population size
	3. T = the time since the population split
	4. m = the symmetric migration rate
	Overview: the ancestral population splits into population 1 with size nu1 and population 2 with size nu2 T generations from the present.
		After the split, there is symmetric migration at rate m.

	[NOTE] all parameters are given relative to the reference size. The ancestral reference size must be calculated from theta
	"""
	# Initialize the parameter values
	nu1, nu2, T, m = params
	# Create an intital 1D AFS that represents our equilibrium ancestral population (use combined sample sizes)
	sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
	# Turn it into a FS object
	fs = moments.Spectrum(sts)
	# Create a 2D FS
	fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
	# After the splitting T generations ago, the populations have constant size and experience symmetric migration
	fs.integrate([nu1, nu2], T, m=np.array([[0, m], [m, 0]]))
	# Return the frequency spectrum obtained from the model
	return fs

# Define the demographic function for a population split followed by exponential change in population 1 and no change in population 2
def split_exp_const_mig(params, ns):
	"""
	Parameter values:
	1. nu1B = the ratio of population 1's bottleneck size to the ancestral population size
	2. nu1C = the ratio of population 1's current size to the ancestral population size
	3. nu2C = the ratio of population 2's current size to the ancestral population size
	4. T = the time since the population split
	5. m = the symmetric migration rate
	Overview: the ancestral population splits into population 1 with size nu1B and population 2 with size nu2C T generations from the present.
		After the split, the population 1 undergoes exponential change to nu1C and population stays constant at nu2C with symmetric migration between them.

	[NOTE] all parameters are given relative to the reference size. The ancestral reference size must be calculated from theta
	"""
	# Initialize the parameter values
	nu1B, nu1C, nu2C, T, m = params
	# Create an intital 1D AFS that represents our equilibrium ancestral population (use combined sample sizes)
	sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
	# Turn it into a FS object
	fs = moments.Spectrum(sts)
	# Create a 2D FS
	fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
	# Create a function that describes exponential growth in pop 1 from size nu1B to nu1C over T generations
	nu1_func = lambda t: nu1B * (nu1C/nu1B) ** (t/T)
	# Make a general nu_func
	nu_func = lambda t: [nu1_func(t), nu2C]
	# After the splitting T generations ago, pop 1 changes exponentially according to function while pop 2 stays the same. And there is symmetric migration between them
	fs.integrate(nu_func, T, m=np.array([[0, m], [m, 0]]))
	# Return the frequency spectrum obtained from the model
	return fs

#-------------------------
# RUN THE MAIN FUNCTION
#-------------------------

if __name__ == '__main__':
	main()




