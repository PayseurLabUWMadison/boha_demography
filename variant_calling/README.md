# Sequence data processing and variant discovery
This Markdown file contains the steps used to process raw sequencing reads, align them to the *P. leucopus* reference genome, and call/filter variants with GATK. 

## Software information
All of the software used is either available in the packages directory or is contained within the `bioinfo` conda environment, which can be constructed from the `bioinfo.yml` file present in the packages directory.

Below are details about each individual software package:
- BBduk (from BBMap version 38.96): available in packages directory
- Trimmomatic (version 0.39): available in packages directory
- BWA-MEM (version 0.7.17-r1188): contained in bioinfo environment
- QualiMap (version 2.2.1): available in the packages directory 
- Picard MarkDuplicates (version 2.25.0 from GATK version 4.2.0.0): available in packages directory
- GATK (version 4.2.0.0): available in packages directory
- FastQC (version 0.11.9): available in packages directory
- Samtools (version 1.7): contained in bioinfo environment

## Methodological overview
What follows are recommended steps that should be taken in processing short-read sequencing data from raw FASTQs to a high confidence variant callset. Much of the information below is based on recommendations from the 2019 version of the [Biostars Handbook](https://www.biostarhandbook.com/) (BSHB), from [Pfeifer (2017)](https://www.nature.com/articles/hdy2016102.pdf), or from GATK's recommendations (with relevant tutorials/guides linked). These steps can broadly involve read quality control, alignment, and variant calling.

**(1) Sequence adapter trimming ([BBduk](https://sourceforge.net/projects/bbmap/)):** this step involves removing Illumina adapter sequences from the reads. These adapter sequences are unique to the sequencing technology used and are required for the DNA fragments to hybridize to the flowcell during short-read sequencing.

**(2) Read quality trimming ([Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)):** this step chews back the ends of reads where the reliability of the base calls decreases. BSHB recommends trying the adapter trimming and quality trimming steps in different orders to see what works best. 

>**Visualization:** this should be performed before and after each QC step to see how the data quality is altered. FASTQ data can be easily visualized with [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/). FastQC is useful for looking at both read quality and adapter contamination. For multiple sequences, it's useful to use the [MultiQC](https://multiqc.info/) tool to aggregate output from FastQC.

>**Removing duplicates:** duplicates refer to reads that were derived from the same exact molecule of DNA. This can happen for a variety of reasons– either PCR amplification bias (during library preparation, the same DNA molecule becomes overrepresented in the pool of PCR products) or due to random chance that two or more PCR products originating from the same fragment happen to hybridize to the flowcell. In either case, this overrepresentation can cause a false sense of reliability in the variant calls at a particular site (i.e., pseudoreplication) or they could cause the propagation of sequencing errors. We actually want to wait to remove duplicates until we've already aligned the reads because if we remove duplicates on the basis of sequence identity and not alignment identify, we could favor erroneous reads.

**(3) Read alignment ([BWA-MEM](http://bio-bwa.sourceforge.net/bwa.shtml)):** this step involves aligning the reads contained within the FASTQ files to the reference genome. BWA-MEM has to first build an index from the reference genome. Then it can be run in paired-end mode to generate the alignments.

>**Visualization:** [QualiMap](http://qualimap.conesalab.org/) is a useful tool for visualizing the quality of alignments, especially when it comes to metric like coverage. 

**(4) Marking duplicates ([Picard MarkDuplicates](https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard-)):** this step involves marking the duplicate reads (as mentioned earlier) that are observed in the alignments. MarkDuplicates will produce new SAM/BAM files with the duplicate reads "tagged". This should be followed by the REMOVE_DUPLICATES option to discard these reads. After the reads have been aligned and the duplicates have been marked, the SAM files should be sorted and compressed into BAM files using samtools. 

>**Visualization:** it's sometimes useful to look at statistics like the percentage of duplicates, read depth, percent unmapped, etc. in the uncompressed SAM files with [samtools](http://www.htslib.org/doc/samtools.html). 

**(5) Variant identification ([GATK HaplotypeCaller](https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller)):** this step involves computing an individual's genotype likelihoods for all possible genotypes at each site from the available read data. HaplotypeCaller's approach is to locally reassemble the haplotypes an individual carries at a locus from their reads, then computing all possible genotype likelihoods. More information about HaplotypeCaller's model can be found at [this](https://gatk.broadinstitute.org/hc/en-us/articles/360035531412) tutorial or at the HaplotypeCaller [preprint](https://www.biorxiv.org/content/10.1101/201178v2.full.pdf). This step results in the creation of a GVCF for each individual, which contains the genotype likelihoods calculated for the entire genome.

**(6) Joint calling ([GATK GentoypeGVCFs](https://gatk.broadinstitute.org/hc/en-us/articles/360037057852-GenotypeGVCFs)):** this step involves combining the genotype likelihoods calculated by HaplotypeCaller for each individual together across all individuals in a populaiton cohort to call variants. By aggregating information across all individuals, this joint-calling approach should [improve the accuracy of the calls](https://gatk.broadinstitute.org/hc/en-us/articles/360035890431-The-logic-of-joint-calling-for-germline-short-variants). In addition to identifying these variants in the sample, this step also generates the QUAL score, which measures our confidence that the site is indeed variable in the sample. More information about the math behind QUAL score calculation can be found [here](https://github.com/broadinstitute/gatk-docs/blob/master/gatk3-methods-and-algorithms/Math_notes:_Understanding_the_QUAL_score_and_its_limitations.md). 

**(7) Variant filtering ([GATK VariantFiltration](https://gatk.broadinstitute.org/hc/en-us/articles/360037434691-VariantFiltration)):** since GATK's approach is meant to be highly sensitive to finding variants, this step is necessary to remove potential false-positive variant calls from the dataset. Although more sophisticated machine learning-based approaches like VQSR are favored over "hard-filtering" variants, these approaches require training data that does not exist for P. leucopus. Because of this, we drew heavily from [this](https://gatk.broadinstitute.org/hc/en-us/articles/360035531112--How-to-Filter-variants-either-with-VQSR-or-by-hard-filtering) tutorial to decide which filtering thresholds to use.

## Pipeline
While the above section outlines the major steps, this section details the exact commands and parameters used. What follows can be divided into three subsections: Read Quality Control, Alignment, and Variant Calling. The following wildcards were used in constructing file names for the example commands:
- Sample name: `{sample}`
- Population name: `{pop}`
- Miscellaneous numbers: `{X}`
- Directory: `{dir}`
- Chromosome: `{chrom}`

#### Read quality control
First, we generated pre-quality control reports with FastQC. This was done for each set of paired-end sequencing files (denoted R1 or R2). At this stage, each sample's paired-end data is represented by multiple sequencing lanes (e.g., L001, L002, etc.) that will depend on the requested depth of coverage.
```
variant_calling/packages/FastQC/fastqc {sample}_S{X}_L00{X}_R1_00{X}.fastq.gz {sample}_S{X}_L00{X}_R2_00{X}.fastq.gz --outdir={dir}
```

Then, we trimmed adapters in all samples with BBduk. This was done using the paired-end setting of BBduk which takes in the forward and reverse read files and creates a trimmed output file for each. 
```
variant_calling/packages/bbmap/bbduk.sh -Xmx1g in1={sample}_S{X}_L00{X}_R1_00{X}.fastq.gz in2={sample}_S{X}_L00{X}_R2_00{X}.fastq.gz out1={sample}_L00{X}_adaptrim_1P.fastq.gz out2={sample}_L00{X}_adaptrim_2P.fastq.gz  ref=variant_calling/reference_files/adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo
```
Afterwards, we created post adapter trimming quality reports with FastQC.
```
variant_calling/packages/FastQC/fastqc {sample}_L00{X}_adaptrim_1P.fastq.gz {sample}_L00{X}_adaptrim_2P.fastq.gz  --outdir={dir}
```

After removing adapters, we trimmed away the low quality ends of reads with Trimmomatic. 
```
java -jar variant_calling/packages/Trimmomatic-0.39/trimmomatic-0.39.jar PE {sample}_L00{X}_adaptrim_1P.fastq.gz {sample}_L00{X}_adaptrim_2P.fastq.gz  -baseout {sample}_L00{X}_adaptrim_qualtrim.fastq.gz SLIDINGWINDOW:4:15 MINLEN:36
```
Then we created post adapter, post quality trimming reports with FastQC.
```
variant_calling/packages/FastQC/fastqc {sample}_L00{X}_adaptrim_qualtrim_1P.fastq.gz {sample}_L00{X}_adaptrim_qualtrim_2P.fastq.gz  --outdir={dir}
```
>Note: This will not improve overall read quality if there are low quality bases in the middle of the read.

#### Alignment
We used the BWA-MEM algorithm to map paired-end reads to the *P. leucopus* reference genome [UCI_PerLeu_2.1](https://www.ncbi.nlm.nih.gov/assembly/GCF_004664715.2). This step was performed separately for each lane's worth of foward and reverse reads.
```
bash variant_calling/scripts/read_group.sh variant_calling/reference_files/GCF_004664715.2_UCI_PerLeu_2.1_genomic.fna.gz {sample}_L00{X}_adaptrim_qualtrim_1P.fastq.gz {sample}_L00{X}_adaptrim_qualtrim_1P.fastq.gz {sample} {sample}_L00{X}_sorted.bam
```
>Note: The `read_group.sh` script is required to properly format the RG tags in the output BAM file. The BWA-MEM algorithm is run from inside this bash script. The BWA command within this script has the following format: `bwa mem -t 8 -R $(echo "@RG\tID:${ID}\tPL:${PL}\tPU:${PU}\tLB:${LB}\tSM:${SM}") $REF $R1 $R2 | samtools sort -o $OUT`.

Then we calculated some quality statistics on each partial alignment file using samtools.
```
samtools flagstat {sample}_L00{X}_sorted.bam > {sample}_L00{X}_sorted.flagstat.txt
```

Then, for each sample, we merged together the lane-separated alignments such that each sample is subsequently represented by a single BAM file. 
```
samtools merge {sample}_sorted_merged.bam {sample}_L00{X}_sorted.bam {sample}_L00{Y}_sorted.bam {sample}_L00{Z}_sorted.bam
```
>Note: The number of additional files following the output file, `{sample}_sorted_merged.bam`, will depend on the number of lanes each individual was sequenced on.

Finally, duplicate reads present in the alignments can be marked for removal with Picard MarkDuplicates. In addition to creating new BAM files with duplicate reads tagged, this also creates a DuplicationMetrics file for each sample.
```
variant_calling/packages/gatk-4.2.0.0/gatk MarkDuplicates --TMP_DIR {dir} -I {sample}_sorted_merged.bam -O {sample}_sorted_merged_rmdups.bam -M {sample}_sorted_merged_DuplicationMetrics.txt
```

#### Variant calling
Before running GATK on the alignments, the BAM files need to be indexed with samtools. At this stage, each individual is represented by a single BAM file.
```
samtools index {sample}_sorted_merged_rmdups.bam
```

To compute genotype likelihoods at all sites and create GVCFs for each individual, we ran HaplotypeCaller on the GVCF Reference Confidence Model mode. This can be done on a per-individual basis, but we found it necessary to additionally split the input BAM for each individual into separate chromosomes, then make the GVCFs on a per-individual, per-chromosome basis.
```
variant_calling/packages/gatk-4.2.0.0/gatk HaplotypeCaller --tmp-dir {dir} -R variant_calling/reference_files/GCF_004664715.2_UCI_PerLeu_2.1_genomic.fna.gz -I {sample}_sorted_merged_rmdups.bam -L {chrom} -ERC GVCF -O {sample}_chr_{chrom}_sorted_merged_rmdups.g.vcf.gz
```
This produces a single GVCF for each individual, for each of their chromosomes.

In order to joint-call variants across individuals from a given population, the per-sample GVCFs need to be combined into a single GVCF that contains all individuals in the population. We did this on a per-chromosome basis for each population using CombineGVCFs.
```
variant_calling/packages/gatk-4.2.0.0/gatk CombineGVCFs --tmp-dir {dir} -R variant_calling/reference_files/GCF_004664715.2_UCI_PerLeu_2.1_genomic.fna.gz -V {sampleX}_chr_{chrom}_sorted_merged_rmdups.g.vcf.gz -V {sampleY}_chr_{chrom}_sorted_merged_rmdups.g.vcf.gz -V {sampleZ}_chr_{chrom}_sorted_merged_rmdups.g.vcf.gz -L {chrom} -O {pop}_chr_{chrom}_sorted_merged_rmdups.g.vcf.gz
```
>Note: The number of input GVCFs specified with each `-V` flag will depend on the number of samples in that population. 
The output of this command is a single combined GVCFs for each chromosome, for each population.

We then performed joint-calling of variants within each population by running GenotypeGVCFs on a per-chromosome basis.
```
variant_calling/packages/gatk-4.2.0.0/gatk GenotypeGVCFs --tmp-dir {dir} -R variant_calling/reference_files/GCF_004664715.2_UCI_PerLeu_2.1_genomic.fna.gz -V {pop}_chr_{chrom}_sorted_merged_rmdups.g.vcf.gz -L {chrom} -O {pop}_chr_{chrom}_sorted_merged_rmdups.vcf.gz
```

Then we merged these per-chromosome unfiltered VCFs into a single file using MergeVCFs.
```
variant_calling/packages/gatk-4.2.0.0/gatk MergeVcfs I={pop}.list O={pop}_full_genome_sorted_merged_rmdups.vcf.gz
```
Here, the `{pop}.list` is a list of all the per-chromosome VCFs for that population, with one file per line.

To obtain a high-confidence variant callset from these unfiltered VCFs, we applied the filters using VariantFiltration. 
```
variant_calling/packages/gatk-4.2.0.0/gatk VariantFiltration -V {pop}_full_genome_sorted_merged_rmdups.vcf.gz -filter "QD < 5.0" --filter-name "QD5" -filter "SOR > 3.0" --filter-name "SOR3" -filter "FS > 60.0" --filter-name "FS60" -filter "MQ < 40.0" --filter-name "MQ40" -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" -O filtered_{pop}_full_genome_sorted_merged_rmdups.vcf.gz
```
>Note: This only tags the records in the FILTER field of the VCF file as "PASS" if they passed all of the filters. This does *not* remove variants by itself, SelectVariants is needed to perform this task.

We then used SelectVariants to restrict the variants contained within each VCF file to only those that passed all of the above filters and were SNP-type variants (as opposed to indels).
```
variant_calling/packages/gatk-4.2.0.0/gatk SelectVariants -V filtered_{pop}_full_genome_sorted_merged_rmdups.vcf.gz --select-type-to-include SNP --exclude-filtered true -O snpOnly_filteredPASS_{pop}_full_genome_sorted_merged_rmdups.vcf.gz
```

We further refined the variants contained in the callsets by also requiring that they be biallelic (as opposed to multiallelic).
```
variant_calling/packages/gatk-4.2.0.0/gatk SelectVariants -V snpOnly_filteredPASS_{pop}_full_genome_sorted_merged_rmdups.vcf.gz --restrict-alleles-to BIALLELIC -O biallelic_snpOnly_filteredPASS_{pop}_full_genome_sorted_merged_rmdups.vcf.gz
```

#### Merging variant callsets
Recall that VCFs for each population will contain sites that are polymorphic only within the observed sample. This complicates SFS-based multipopulation demographic inference because it limits our ability to say whether a site observed as polymorphic in one sample is monomorphic vs uncalled in another sample. The steps below describe how we "rescued" the genotype calls at these "private" variants in order to properly merge the callsets across all three populations.

First we need to identify sites that are unique to each population. This can be accomplished with bcftools `isec`. For each focal population, we start by identifying records present in each of the *other* two populations, but *not* in our focal population.

This command prints records present in pop1 but *not* our focal population:
```
bcftools isec -C -Oz -o {focal}_unique_to_{pop1}.txt biallelic_snpOnly_filteredPASS_{pop1}_full_genome_sorted_merged_rmdups.vcf.gz biallelic_snpOnly_filteredPASS_{focal}_full_genome_sorted_merged_rmdups.vcf.gz
```

This command prints records present in pop2 but *not* our focal population:
```
bcftools isec -C -Oz -o {focal}_unique_to_{pop2}.txt biallelic_snpOnly_filteredPASS_{pop2}_full_genome_sorted_merged_rmdups.vcf.gz biallelic_snpOnly_filteredPASS_{focal}_full_genome_sorted_merged_rmdups.vcf.gz
```

We then merge these lists of private SNPs together:
```
cat {focal}_unique_to_{pop1}.txt {focal}_unique_to_{pop2}.txt | cut -f1,2 | sort -k1,1 -k2,2n | uniq > {focal}_to_force_call.txt
```
And use this combined list to create a [GATK-compatible interval file](https://gatk.broadinstitute.org/hc/en-us/articles/360035531852-Intervals-and-interval-lists):
```
awk '{print $1":"$2"-"$2}' {focal}_to_force_call.txt > {focal}_to_force_call.list
```

We create this `{focal}_to_force_call.list` for *each* of the three populations.

As with the previous variant calling steps, we found it necessary to perform these steps in a per-chromosme manner. This starts with splitting the `{focal}_to_force_call.list` by chromosome:
```
cat {pop}_to_force_call.list | grep {chrom} > {chrom}_{popl}_to_force_call.list
```

We can then pass these lists of variants that are "missing" from our focal population callset to GATK and use the `-all-sites` option to output records for these invariant sites:
```
variant_calling/packages/gatk-4.2.0.0/gatk GenotypeGVCFs -R variant_calling/reference_files/GCF_004664715.2_UCI_PerLeu_2.1_genomic.fna.gz -V {pop}_chr_{chrom}_sorted_merged_rmdups.g.vcf.gz -L {chrom}_{pop}_to_force_call.list -ip 100 -all-sites -stand-call-conf 0.0 -O {chrom}_{pop}_force_calls.vcf.gz
```

We still want to filter these "rescued" records to ensure that we have confidence in the genotypes being called. However, in contrast to the filtering performed above, we omit any QUAL score-based filters since, by nature of these sites being invariant, they will yield poor QUAL scores:
```
variant_calling/packages/gatk-4.2.0.0/gatk VariantFiltration -V {chrom}_{pop}_force_calls.vcf.gz -filter "SOR > 3.0" --filter-name "SOR3" -filter "FS > 60.0" --filter-name "FS60" -filter "MQ < 40.0" --filter-name "MQ40" -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" -O filtered_{chrom}_{pop}_force_calls.vcf.gz
```

Then we exclude those that failed the above filters and remove and sites that are neither SNPs nor invariant sites (e.g., indels) using SelectVariants:
```
variant_calling/packages/gatk-4.2.0.0/gatk SelectVariants -V filtered_{chrom}_{pop}_force_calls.vcf.gz --select-type-to-include NO_VARIATION --select-type-to-include SNP --exclude-filtered true -L {chrom}_{pop}_to_force_call.list -O snpOnly_filteredPASS_{chrom}_{pop}_force_calls.vcf.gz
```

Once these invariant calls have been "rescued", we can combine them with the original variant callset produced for each population (expanding {chromX}, {chromY}, and {chromZ} to reflect all of the chromosome-split VCFs):
```
variant_calling/packages/gatk-4.2.0.0/gatk MergeVcfs I=biallelic_snpOnly_filteredPASS_{pop}_full_genome_sorted_merged_rmdups.vcf.gz I=snpOnly_filteredPASS_{chromX}_{pop}_force_calls.vcf.gz I=snpOnly_filteredPASS_{chromY}_{pop}_force_calls.vcf.gz I=snpOnly_filteredPASS_{chromZ}_{pop}_force_calls.vcf.gz O=allsites_snpOnly_filteredPASS_{pop}_full_genome_sorted_merged_rmdups.vcf.gz
```

After constructing these monomorphic+polymorphic callsets for each population, we can combine them to create a multi-population callset using bcftools `merge`. 
```
bcftools merge -m all -Oz -o allsites_snpOnly_filteredPASS_boha_full_genome_sorted_merged_rmdups.vcf.gz allsites_snpOnly_filteredPASS_{pop1}_full_genome_sorted_merged_rmdups.vcf.gz allsites_snpOnly_filteredPASS_{pop2}_full_genome_sorted_merged_rmdups.vcf.gz allsites_snpOnly_filteredPASS_{pop3}_full_genome_sorted_merged_rmdups.vcf.gz
``` 

Finally, to remove any "rescued" sites that are monomorphic in *all* populations, we can use GATK SelectVariants to restrict the combined callset to only include sites that are biallelic and polymorphic across the three populations.
```
variant_calling/packages/gatk-4.2.0.0/gatk SelectVariants -V allsites_snpOnly_filteredPASS_boha_full_genome_sorted_merged_rmdups.vcf.gz --select-type-to-include SNP --restrict-alleles-to BIALLELIC -O biallelic_allsites_snpOnly_filteredPASS_boha_full_genome_sorted_merged_rmdups.vcf.gz
```

## Additional filtering procedures
The above steps yield a multi-population callset that has been filtered for variant quality. Additional filtering procedures used to exclude SNPs falling in inaccessible or genic regions are described in the Materials and Methods section. 


