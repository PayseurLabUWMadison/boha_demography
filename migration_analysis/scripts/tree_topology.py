#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# OVERVIEW
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Purpose: analyze the results of a given replicate tree sequence (note that each replicate represents a 100 Mbp chromosome)
# Output: a single data frame that contains, for each tree in the sequence:
# [chrom] [model] [mig_rate] [mig_time] [genomic coordinates] [true/false whether it involves migration] [tree height] [summed branch lengths] [sum internal branch lengths] [sum external branch lengths] [sum of basal branch lengths] [average external branch length] [variance in external branch length]
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Run as: 
# python tree_topology.py --ts [ts file] --model [model name] --mig_rate [mig rate] --mig_time [mig time] --chrom [chrom]
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------

#---------------------------
# IMPORT PACKAGES
#---------------------------
import msprime
import tskit
import pandas as pd
import numpy as np
import argparse
import statistics

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

	# Tree sequence to analyze
	tsf = args.ts
	# Set the demographic model used to simulate the data
	model = args.model
	# Set the migration rate used to simulate the data
	mig_rate = float(args.mig_rate)
	# Set the migration timing (as a fraction of total split time) used to simulate the data
	mig_time = float(args.mig_time)
	# Specifies the chromosome used to simulate the data
	chrom = int(args.chrom)

	# This is the name of the topology file
	topo_out = "Model-" + str(args.model) + "_MigRate-" + str(args.mig_rate) + "_MigTime-" + str(args.mig_time) + "_Chrom-" + str(args.chrom) + ".topology"
	# This is the name of the migration tract file
	mig_out = "Model-" + str(args.model) + "_MigRate-" + str(args.mig_rate) + "_MigTime-" + str(args.mig_time) + "_Chrom-" + str(args.chrom) + ".migrant_tracts"

	print("Done.")

	#----------------------
	# Iterate over trees
	#----------------------

	print("Loading the tree sequence...")

	# Load the tree sequence
	ts = tskit.load(tsf)

	# Create empty lists to store relevant information
	start_ls = []
	stop_ls = []
	mig_bool_ls = []
	mig_hap_ls = []
	height_ls = []
	sum_bl_ls = []
	sum_inbl_ls = []
	sum_exbl_ls = []
	sum_babl_ls = []
	avg_exbl_ls = []
	sd_exbl_ls = []
	root_pop_ls = []

	print("Done.")

	print("Collecting information from the tree sequence...")

	# Create variable to store source pop id
	source_pop = ""

	# Loop through populations
	for i in ts.populations():
		# Print out info
		print("Population " + str(i.id) + " is named " + str(i.metadata['name']))
		# Determine which of the pop ids represents the ghost population
		if (i.metadata['name'] == "Ghost"):
			# Store id of ghost population
			source_pop = i.id
		# Determine if this id represents the ancestral pop
		elif (i.metadata['name'] == "Ancestral"):
			# Store the id
			anc_pop = i.id

	# Create a list to store migrant tract info (as left pos, right pos, time)
	migrating_tracts = []

	# Loop through each entry in the migration table
	for migration in ts.migrations():
		# Determine whether this tract came from the ghost population (i.e., is the "destination" going back in time)
		if migration.dest == source_pop:
			# If so, append the left and right coordinates along with the time
			migrating_tracts.append((migration.left, migration.right, migration.time))

	# Print info aboout migrations
	print("Migration table:")
	print(ts.tables.migrations)

	print("ID of source population for migrant segments is: " + str(source_pop))

	# Print the number of migrations observed in this tree sequence
	print("This tree sequence has " + str(ts.num_migrations) + " total migrations")
	print(str(len(migrating_tracts)) + " of those are from the source population")
	print(migrating_tracts)
	print("This tree sequence has " + str(ts.num_trees) + " total trees")
	
	print("Writing migrant tract output")

	# Write the output of the migration tables
	df = pd.DataFrame(migrating_tracts, columns = ["start", "stop", "time"])
	# Append additional info
	df["model"] = model
	df["mig_rate"] = mig_rate
	df["mig_time"] = mig_time
	df["chrom"] = chrom
	# Write the output to tab delim csv
	df.to_csv(mig_out, sep="\t", index=False)

	print("Done.")

	print("Analyzing tree topologies...")

	# Loop through the list of all trees
	for tree in ts.trees():
		# # Print info about this tree
		# print(f"Tree {tree.index} covers {tree.interval}")

		# Indicator variable for whether this tree falls within migrant segment
		mig_bool = False
		# Counter to keep tract of the number of migrant haplotypes 
		mig_hap = 0

		# Check to see if the interval for this tree falls within any of the migrant segments
		for tract in migrating_tracts:
			# Gather the left pos, right pos, and time from this entry in the full list
			mig_left = tract[0]
			mig_right = tract[1]
			# Check to see if the interval for this tree falls within given segment
			if (tree.interval.left >= mig_left and tree.interval.right <= mig_right):
				# And change indicator variable
				mig_bool = True
				# And add one to the number of migrant haplotypes in the sample at this tree
				mig_hap = mig_hap + 1

		# Once we've checked all the migrant tracts, append the result for this tree to the lists
		mig_bool_ls.append(mig_bool)
		mig_hap_ls.append(mig_hap)

		# Print progress report
		if (tree.index % 1000 == 0):
			print("Analyzed " + str(tree.index) + " tree sequences")

		# Compute the genomic coordinates of this tree
		start_ls.append(tree.interval.left)
		stop_ls.append(tree.interval.right)
		# Compute the sum of branch lengths
		sum_bl_ls.append(tree.total_branch_length)

		# Check to see if there is more than one root for this tree
		if (tree.has_single_root):
			# # Print info
			# print("This tree has one root")
			# Store the node ID of the root
			root = tree.root
			# Figure out what the population id is for the root node
			root_pop_id = tree.population(root)
			# Use the population table to figure out what the name is
			root_pop_name = ts.population(root_pop_id).metadata['name']
			# And append this to the list
			root_pop_ls.append(root_pop_name)
			# Get the topological stats from the function
			height,sum_exbl,avg_exbl,sd_exbl,sum_inbl,sum_babl,exbls = topology(tree, root)
			# Append these values to the lists
			height_ls.append(height)
			sum_exbl_ls.append(sum_exbl)
			avg_exbl_ls.append(avg_exbl)
			sd_exbl_ls.append(sd_exbl)
			sum_inbl_ls.append(sum_inbl)
			sum_babl_ls.append(sum_babl)

		# In case the given tree has multiple roots
		else:
			# # Print info
			# print("This tree has multiple roots")
			# Store all the root IDs
			roots = tree.roots
			# Append filler value for root pop name
			root_pop_ls.append("NA")
			# Create lists to hold values for each root
			node_height_ls = []
			node_sum_inbl_ls = []
			node_sum_exbl_ls = []
			node_sum_babl_ls = []
			node_avg_exbl_ls = []
			node_exbl_ls = []
			# Perform the topological analysis for each root
			for root in roots:
				# Get the topological stats from the function
				height,sum_exbl,avg_exbl,sd_exbl,sum_inbl,sum_babl,exbls = topology(tree, root)
				# Append these values to the lists
				node_height_ls.append(height)
				node_sum_exbl_ls.append(sum_exbl)
				node_avg_exbl_ls.append(avg_exbl)
				node_exbl_ls.append(exbls)
				node_sum_inbl_ls.append(sum_inbl)
				node_sum_babl_ls.append(sum_babl)
			# Sum and re-average branch lengths across nodes
			height_ls.append(max(node_height_ls))
			sum_inbl_ls.append(sum(node_sum_inbl_ls))
			sum_exbl_ls.append(sum(node_sum_exbl_ls))
			sum_babl_ls.append(sum(node_sum_babl_ls))
			avg_exbl_ls.append(statistics.mean(node_avg_exbl_ls))
			sd_exbl_ls.append(statistics.stdev(node_exbl_ls))

		# # Break statement for testing purposes
		# if (tree.index == 5000):
		# 	break

	# Once we're done constructing all the neccesary lists, turn it into a dataframe
	df = pd.DataFrame({"model": model, "mig_rate": mig_rate, "mig_time": mig_time, "chrom": chrom, "start": start_ls, "stop": stop_ls, "mig": mig_bool_ls, "num_mig_haps": mig_hap_ls, "height": height_ls, "total_branch_lengths": sum_bl_ls, "sum_internal_branch_lengths": sum_inbl_ls, "sum_external_branch_lengths": sum_exbl_ls, "sum_basal_branch_lengths": sum_babl_ls, "mean_external_branch_lengths": avg_exbl_ls, "sd_external_branch_lengths": sd_exbl_ls, "mrca_pop": root_pop_ls})
	# print(df)

	print("Writing topology output...")
	df.to_csv(topo_out, sep="\t", index=False)

#---------------------
# DEFINE FUNCTIONS
#---------------------

def topology(tree, root):
	"""
	For the given tree, compute a variety of branch length stats for the specified root
	Need to return the full exbl list to recompute sd_exbl in the main function for the case of multiple roots
	returns: height,sum_exbl,avg_exbl,sd_exbl,sum_inbl,sum_babl,exbls
	"""
	# Store the nodes that are reachable from this single root
	nodes = tree.nodes(root)
	# print(root)
	# Compute the tmrca based on the first node
	height = tree.tmrca(0, root)
	# print(height)

	# Create a list to store leaf nodes
	leaves = []
	# Create a list to store internal nodes
	internals = []
	# Create a list to store the nodes whose parent is the root
	basals = []

	# Loop through the nodes and determine whether they are internal nodes or external "leaf" nodes
	for node in nodes:

		# Check to see whether this node is a leaf or internal and append to corresponding list
		if (tree.is_leaf(node)):
			leaves.append(node)
		else:
			internals.append(node)

		# Additionally check to see if the node's parent is the root (i.e., it is a basal node)
		# Note that the internal branch lengths or external branch lengths could include the basal branch lengths
		if (tree.parent(node)==root):
			basals.append(node)

	# print(leaves)
	# print(internals)
	# print(basals)

	# Using the leaf node list, grab corresponding branch lengths
	exbls = [tree.branch_length(x) for x in leaves]
	# From these external branch lengths, compute the mean, sd, and sum of branch lengths
	sum_exbl = sum(exbls)
	avg_exbl = statistics.mean(exbls)
	sd_exbl = statistics.stdev(exbls)

	# Using the internal node list, grab corresponding branch lengths
	inbls = [tree.branch_length(x) for x in internals]
	# From these internal branch lengths, compute the sum of branch lengths
	sum_inbl = sum(inbls)

	# Using the basal node list, grab the corresponding branch lenths
	babls = [tree.branch_length(x) for x in basals]
	# From these basal branch lenths, compute the sum of branch lengths
	sum_babl = sum(babls)

	return height,sum_exbl,avg_exbl,sd_exbl,sum_inbl,sum_babl,exbls

def parse_args():
	"""
	Parse the command-line arguments
	"""
	# Call to argparse
	parser = argparse.ArgumentParser()

	# Command-line variables

	# Tree sequence input
	parser.add_argument('--ts', required=True)	
	# Determines which demographic model should be used
	parser.add_argument('--model', required=True)
	# Determines what the symmetric migration rate parameter should be
	parser.add_argument('--mig_rate', required=True)
	# Determines when (as a fraction of total split time) migration should begin
	parser.add_argument('--mig_time', required=True)
	# Determines which chromosome is being simulated
	parser.add_argument('--chrom', required=True)

	# Parse and return arguments
	args = parser.parse_args()
	return args

#-------------------------
# RUN THE MAIN FUNCTION
#-------------------------

if __name__ == '__main__':
	main()






