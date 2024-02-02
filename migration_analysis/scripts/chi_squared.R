#!/usr/bin/env Rscript
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# OVERVIEW
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Purpose: to compute the chi-squared test statistic on individual singleton counts for each chromosome
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Run as: 
# Rscript chi_squared.R [singleton file] [popfile] [output prefix]
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Load libraries
library(dplyr)

# Read command line input
args <- commandArgs(trailingOnly=TRUE)
singletons <- args[1]
popfile <- args[2]
output <- args[3]

# Read in the singleton output and assign column headers
sf <- read.table(singletons)
colnames(sf) <- c("chrom","pos","designation","indv","population","params")

# Read in the popfile output and assign column headers
pf <- read.table(popfile)
colnames(pf) <- c("indv","population")

# Compute the total sample sizes observed for each population and append to the singleton df
pops <- unique(sf$population)
for(pop in pops){
  nindv <- nrow(pf[pf$population==pop,])
  sf$tot_indv[sf$population==pop] <- nindv
}

# Function to compute the chi_squared statistic for the given table
get_chi_squared <- function(sub){
  
  # Compute the expected number of individuals from the full table's deg free column
  num_indv <- unique(sub$tot_indv)
  
  # Create a table to hold summed singleton observations for each individual
  tab <- as.data.frame(table(sub$indv))
  colnames(tab) <- c("indv", "obs")
  
  # Check to make sure that there are the correct number of individuals
  indv_count <- length(tab$indv)
  
  # If there aren't the correct number of individuals
  if(indv_count < num_indv){
    # Compute the differece
    diff <- num_indv - indv_count
    # Make placeholder individuals and observations
    placeholders <- seq(1:diff)
    indv <- c(as.character(tab$indv), placeholders)
    obs <- c(as.numeric(tab$obs), rep(0, diff))
    # Construct a new replacement df
    tab <- data.frame(indv, obs)
  }
  
  # Create a table to hold expected singleton observations for each individual
  N <- sum(tab$obs)
  prob <- 1/num_indv
  tab$exp <- N*prob
  
  # Compute the chi_squared test statistic
  chi_squared <- sum(((tab$obs-tab$exp)^2)/tab$exp)
  return(chi_squared)
}

# Perform the chi squared analysis for each chromosome
chrom_analysis <- sf %>% group_by(params, chrom, population, tot_indv) %>% do(data.frame(chi_squared=get_chi_squared(.)))
write.table(chrom_analysis, file=paste(output, ".chi_squared", sep=""), quote=FALSE, append=FALSE, row.names=FALSE)

