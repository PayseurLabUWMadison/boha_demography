# Migration simulations and singleton analysis
This Markdown file details the steps used to augment the best-fit demographic models with various scenarios of migration between focal populations and a single, unsampled "ghost" population. Following these simulations, we compared the genomic and sample distribution of singletons between the empirical island datasets and simulated datasets that either included or excluded migration.

## Software information
All of the software used is contained within the `singletons` conda environment, which can be constructed from the singletons.yml file present in the packages directory.

Below are details about individual software packages:
- msprime (version 1.2.0)
- vcftools (version 0.1.16)

## Migration simulations
In order to characterize the effect of migration on the spatial (i.e., along the genome), and sample (i.e., among individuals) distribution of singleton variants, we can conduct chromosome-scale simulations using the inferred best-fit demographic parameters as a starting point and adding various forms of migration. The `sim_migration.py` script in the scripts directory acheives this by taking in a user-specified migration rate and migration duration, then simulating chromosome-sized genomic elements under the given demographic model. To increase computational efficiency, we ran this simulation script in a parallel fashion using UW-Madison's high-throughput computing system in order to obtain a genome-scale dataset of 25 100 Mbp "chromosomes" for each model/migration regime. Details about the migration parameters used can be found in the Materials and Methods section of the manuscript. The output of this script is a VCF-like file for a single genomic element under the specified demographic model. These individual chromosome VCFs can either be analyzed separately or concatenated within models to produce a genome-scale VCF dataset.

## Singleton analysis
In order to determine the location of singleton variants along the genome and the individual that each singleton arose in, we can use the vcftools `--singleton` function. The `analyze_singletons.py` script serves as a wrapper script for the vcftools function (which is run inside of the `vcftools_singletons.sh` that is also available in the scripts directory). The output of this script is a tab-delimited data frame that contains `[chrom] [pos] [designation] [indv] [population] [params]` where the designation refers to either true singletons (S) or private doubletons (D) and the indv refers to the individual in which the variant was observed. Since both the simulated and empirical VCF files on which this script is run contain multiple populations, it's neccessary to pass a population manifest to this script with the `--popfile` command that contains each individual's associated population in the format `[indv] [population]`. 

In addition to allowing us to plot the location of singleton variants along the genome, this output also enables us to test whether singleton variants are uniformly distributed across individuals (i.e., does each individual contribute an equal proportion of the total singletons observed in the sample?). The `chi_squared.R` script in the scripts directory computes the chi-squared test statistic for the `analyze_singletons.py` output for each simulated or observed chromosome to quantify the deviation of individual singleton counts from the null hypothesis that singletons are contributed equally among individuals. This Rscript takes in the output of the `analyze_singletons.py` script and the associated population manifest and outputs the chi-squared test statistic computed at the chromosome level. 











