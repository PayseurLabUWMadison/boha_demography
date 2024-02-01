# Demographic inference and predictive simulations
This Markdown file details the steps used to fit two-population demographic models to the jSFS, perform model comparisons, estimate parameter uncertainties, simulate sequence data under the best-fit demographic models, and compute windowed summary statistics from both simulated and empirical data.

## Software information
All of the software used is contained within the `demog` conda environment, which can be constructed from the [demog.yml](https://github.com/PayseurLabUWMadison/boha_demography/blob/main/demographic_inference/packages/demog.yml) file present in the [packages](https://github.com/PayseurLabUWMadison/boha_demography/tree/main/demographic_inference/packages) directory.

Below are details about individual software packages:
- moments (version 1.1.15)
- msprime (version 1.2.0)
- scikit-allel (version 1.3.6)
- dadi (version 2.2.0)

## Demographic inference

#### Model specification, fitting, and optimization
The first step in demographic model fitting is to generate a folded jSFS for each pair of populations being analyzed. This is acheived using [`create_fs.py`](https://github.com/PayseurLabUWMadison/boha_demography/blob/main/demographic_inference/scripts/create_fs.py) in the [scripts](https://github.com/PayseurLabUWMadison/boha_demography/tree/main/demographic_inference/scripts) directory. This script takes in a filtered, multi-population VCF file and outputs a `.fs` object that can be read into downstream moments scripts. 

After creating the frequency spectrum object, we can fit pre-specified demographic models and find their maximum likelihood parameter estimates. This is acheived using [`fit_models.py`](https://github.com/PayseurLabUWMadison/boha_demography/blob/main/demographic_inference/scripts/fit_models.py) in the [scripts](https://github.com/PayseurLabUWMadison/boha_demography/tree/main/demographic_inference/scripts) directory. This script takes in the two-population `.fs` file, a name to ID the population comparison, the name of the model being fit, a series of optimization parameters, and the effective sequence length and mutation rate used to convert parameter estimates. Based on the optimization parameters provided, the script will then identify the maximum likelihood parameter estimates for the specified model and output them as a tab-delimited file. To ensure that the starting parameter space is thoroughly searched, we ran multiple independent instances of this script in parallel using UW Madison's high-throughput computing system.

#### Model comparison
Once we obtain the maximum likelihood parameter estimates for each model, we can inspect the fit of the model to both the joint and marginal SFS of each population by computing the residual differences. This is acheived using [`model_eval.py`](https://github.com/PayseurLabUWMadison/boha_demography/blob/main/demographic_inference/scripts/model_eval.py) in the [scripts](https://github.com/PayseurLabUWMadison/boha_demography/tree/main/demographic_inference/scripts) directory. This script takes in the best-fit (un-converted) parameter estimates, the model they correspond to, and the `.fs` of the given population comparison then computes and outputs the associated 1D and 2D residual plots. 

For nested models that showed good parameter convergence, we used adjusted likelihood ratio tests to select between contrasting models. These comparisons are hard-coded in the [`model_lrt.py`](https://github.com/PayseurLabUWMadison/boha_demography/blob/main/demographic_inference/scripts/model_lrt.py) script included in the [scripts](https://github.com/PayseurLabUWMadison/boha_demography/tree/main/demographic_inference/scripts) directory. 

#### Estimating parameter uncertainties
Once we've identified best-fit models for each two-population comparison, we can employ the GIM methods developed by [Coffman et al. (2016)](https://academic.oup.com/mbe/article/33/2/591/2579696) to estimate parameter uncertainties. This analysis is coded in the [`est_uncert.py`](https://github.com/PayseurLabUWMadison/boha_demography/blob/main/demographic_inference/scripts/est_uncert.py) script in the [scripts](https://github.com/PayseurLabUWMadison/boha_demography/tree/main/demographic_inference/scripts) directory. This script takes in the best-fitting parameter estimates for the simple split model, creates bootstrap frequency spectra from the input VCF, and uses these bootstrap replicates to compute the GIM and obtain standard deviations for each parameter. 

## Predictive simulations

#### Simulating the best-fit demographic model
In order to see how well our inferred model recovers other features of neutral variation, we can simulate sequence data under the best-fit model for each pair of populations. This is acheived using [`sim_models.py`](https://github.com/PayseurLabUWMadison/boha_demography/blob/main/demographic_inference/scripts/sim_models.py) in the [scripts](https://github.com/PayseurLabUWMadison/boha_demography/tree/main/demographic_inference/scripts) directory. Each of the best-fit demographic models are specified in this script using msprime, with maximum likelihood estimates for relevant parameters hard-coded. Based on user-specified simulation parameters, this script then generates VCF-like data for the given populations/model.

#### Computing summary statistics from simulated and observed data
After we've simulated sequences under the best-fit models, we can compare how summary statistics produced by these inferred histories compare to estimates from the observed data. This is acheived separately for the simulated and observed data using [`calc_stats_simdata.py`](https://github.com/PayseurLabUWMadison/boha_demography/blob/main/demographic_inference/scripts/calc_stats_simdata.py) and [`calc_stats_obsdata.py`](https://github.com/PayseurLabUWMadison/boha_demography/blob/main/demographic_inference/scripts/calc_stats_obsdata.py), respectively. These scripts take in the corresponding VCF file and output pi, Tajima's D, and Fst computed in 5kbp non-overlapping windows. 

The difference between these two scripts comes from the use of a mask file for the observed sequence data. Because of the extensive filtering performed on the empirical sequence data to exclude low-quality SNPs and those falling in "inaccessible" or genic regions, each window will contain a variable number of "accessible" sites. For per-site summary statistics like pi, this heterogeneity has to be accounted for. The mask file indicates the accessibility status for all positions in a given chromosome/contig and is then passed to the `is_accessible` argument of relevant scikit-allel functions (such as `allel.windowed_diversity()`). 

In order for this mask file to be used properly by the `calc_stats_obsdata.py` script, it needs to be in FASTA file format (with separate headers for each chromosome/contig) and inaccessible sites should be replaced with the character `N`. Given a BED file that contains the sites/intervals excluded after the quality, callability, and neutrality-based filtering described in the Materials and Methods section, we can use bedtools `maskfasta` to create a masked-version of the reference genome:
```
bedtools maskfasta -fi [reference FASTA] -bed [BED file of excluded sites/intervals] -fo [masked FASTA]
```
