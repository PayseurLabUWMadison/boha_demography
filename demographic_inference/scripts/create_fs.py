#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# OVERVIEW
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Purpose: the purpose of this script is to generate a moments FS object from an input VCF file
# Output: a moments FS object
# Special inputs: --popfile is the corresponding sample manifest in the format [individual] [pop]
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Run as: 
# python create_fs.py --vcf [filtered VCF file] --fspre [prefix for FS file] --popfile [corresponding sample manifest] --dimension [1, 2, or 3]
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
	# This is the name of the output FS object
	fs_name = args.fspre + ".fs"
	# Variable for holding input popfile
	popfile = args.popfile
	# Create a variable to hold the number of dimensions
	dim = int(args.dimensions)

	print("Done.")

	#-----------------------------------------
	# Constructing frequency spectrum objects
	#-----------------------------------------

	# Run function to create FS from the input VCF
	print("Creating moments frequency spectrum object...")
	create_fs(vcf, fs_name, popfile, dim)
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

	# Determines which VCF to use
	parser.add_argument('--vcf', required=True)
	# Prefix of FS files
	parser.add_argument('--fspre', required=True)
	# Determines which popfile to use
	parser.add_argument('--popfile', required=True)
	# This is the number of dimensions of the input FS
	parser.add_argument("--dimensions", required=True)

	# Parse and return arguments
	args = parser.parse_args()
	return args

def create_fs(vcf, fs_name, popfile, dim):
	"""
	This function takes in a VCF dataset and creates a moments-fs object, outputting this to a loadable .fs file and plotting in a PDF

	vcf: this is the input data
	fs_name: this is the name given to the FS fileset
	popfile: this is the corresponding sample manifest in the format [individual] [pop]

	return: NULL
	"""
	# Read in the input vcf and construct moments dictionary object
	dd = moments.Misc.make_data_dict_vcf(vcf, popfile)

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

	# Create a frequency spectrum for the number of populations specified in the popfile
	# Don't project down samples (i.e. not considering SNPs with missing data)
	# Don't polarize SNPs (naive of ancestral state, means that the FS will be folded)
	fs = moments.Spectrum.from_data_dict(dd, pop_ls, projections = proj_ls, polarized = False)

	# Write the frequency spectrum to a file
	fs.to_file(fs_name)

	# For 1 dimensional analysis
	if (dim==1):
		# Open the figure
		pylab.figure(figsize=(7,6))
		# Plot the folded SFS and save as a PDF
		plot_1d_fs(fs, pop_ls[0])
		# Save the figure
		pyplot.show()
		pyplot.savefig((fs_name + ".pdf"))

	# For 2 dimensional analysis
	elif (dim==2):
		# Open the figure
		pylab.figure(figsize=(7,6))
		# Plot the folded SFS and save as a PDF
		moments.Plotting.plot_single_2d_sfs(fs)
		# Save the figure
		pyplot.show()
		pyplot.savefig((fs_name + ".pdf"))

def plot_1d_fs(fs, name, fig_num=None, show=True):
	"""
	Plot a 1-dimensional frequency spectrum.

	fs: 1-dimensional Spectrum
	fig_num: Clear and use figure fig_num for display. If None, an new figure
			 window is created.
	show: If True, execute pylab.show command to make sure plot displays.

	Note that all the plotting is done with pylab. To see additional pylab
	methods: "import pylab; help(pylab)". Pylab's many functions are documented
	at http://matplotlib.sourceforge.net/contents.html
	"""
	if fig_num is None:
		fig = pylab.gcf()
	else:
		fig = pylab.figure(fig_num, figsize=(7,7))
	fig.clear()

	ax = fig.add_subplot(1,1,1)
	ax.semilogy(fs, '-ob')

	ax.set_xlim(0, fs.sample_sizes[0]-(fs.sample_sizes[0]/2-1))
	fig.suptitle(str(name))
	ax.set_xlabel('Folded SFS bin')
	ax.set_ylabel('Number of SNPs')
	if show:
		fig.show()

#-------------------------
# RUN THE MAIN FUNCTION
#-------------------------

if __name__ == '__main__':
	main()




