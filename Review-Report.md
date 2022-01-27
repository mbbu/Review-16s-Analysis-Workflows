# 16S-Report-Prep
## Introduction
16S rRNA is a gene that encodes the RNA component of the small subunit(30S subunit) of ribosomes in bacteria and archaea. 16S is a sedimentation coefficient([Dependency Map of Proteins in the Small Ribosomal Subunit](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.0020010)). This is an essential gene required for initiating protein synthesis and the stabilizing correct codon-anticodon pairing in the A site of the ribosome during mRNA translation([The distribution, diversity, and importance of 16S rRNA gene introns in the order Thermoproteales](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4496867/)). It is called "the molecular fossil" of bacteria because of being highly conserved and specific. This makes it the most widely used gene marker for genus and species identification, as well as taxonomic significance([16S/18S/ITS Amplicon Sequencing](https://www.cd-genomics.com/16S-18S-ITS-Amplicon-Sequencing.html)). The gene is about 1500bp and is composed of both conserved regions and variable regions. The conserved region is shared while the variable regions have differences among different bacteria and therefore providing information on the specificity of the genus and the species([16S rRNA, One of the Most Important rRNAs](https://www.cd-genomics.com/blog/16s-rrna-one-of-the-most-important-rrnas/)). 

The gene is used in microbiome analysis. Analysis pipelines have been developed and improved over the years and include QIIME and DADA2 pipelines. Our internship project required us to review some of the existing pipelines and come up with a conclusion on the best ones. After testing them with different datasets, we concluded that the nf-core/ampliseq pipeline was the best and also made a few changes in the MBBU 16S-Accreditation pipelines (Both Dada2 and Qiime2 pipelines). Those that were challenging to run were dropped early in the review process.

## Objectives

* To review existing microbiome workflows, identify great ones and extend the workflows where there are gaps, especially to make them useful in insect and pathogen data.

## Methods

1. Testing workflows
* We tested the following pipelines according to the functionalities indicated in their documentations:
  - https://github.com/nf-core/ampliseq 
  - https://h3abionet.github.io/H3ABionet-SOPs/16s-rRNA-1-0.html
  - https://github.com/h3abionet/TADA
  - https://github.com/mbbu/16S_Accreditation
  - https://github.com/h3abionet/h3abionet16S
  - Yosef's DADA2 pipeline
* We used different datasets in running the pieplines

2. Identifying gaps
* While testing the pipelines, we identified gaps using the following criteria:
  - How easy are they to set up and use? Do they provide accessible documentation and tutorials?
  - Are they fast and easily scalable based on available compute resources?
  - Can they scale to the cloud?
  - Can they be used on a variety of data, including insects and pathogen microbiome
  - Are they implemented in the latest specifications and versions of the tools? For example, Whether the pipeline implements Nextflow DSL2 syntax and docker or singularity containers
  - Are they well and regularly maintained? When were they updated last?
* We were also able to find gaps from the errors we came across

3. Extending workflows
* We worked around the errors we got and extensively tested the final extended workflows using different sets of data

## Results
* [The new MBBU/16S-Accreditation DADA2 pipeline](https://github.com/mbbu/Reviewing-16s-Analysis-Workflows/tree/main/New_MBBU_Dada2_pipeline)

## Summary of work done
| WEEK | ACTIVITY |
|----- | -------- |
| WEEK 1 | Obtaining test datasets 
|        | Assessing workflow performance 
| WEEK 2 | Running and troubleshooting MBBU/16S-Accreditation, H3ABioNet-SOPs, H3ABioNet-TADA, H3ABioNet-16S, nf-core/ampliseq
| WEEK 3 | Testing the nf-core/ ampliseq using stingless-bee data
|        | Running Yosef's DADA2 pipeline and MBBU/16S-Accreditation DADA2 pipeline
| WEEK 4 | Obtaining more test data
|        | Testing of the nf-core pipeline using different datasets(ITS data and 18S data)
|        | Solving week 3 errors and running Yosef's DADA2 pipeline, MBBU/16S-Accreditation DADA2 pipeline
| WEEK 5 | Solving Week 4 errors and running Yosef's DADA2 pipeline, MBBU/16S-Accreditation DADA2 pipeline
|        | Testing of the MBBU/16S-Accreditation DADA2 pipeline using different datasets(stingless bee, dog, and nf-core/ampliseq data)
|        | Testing of the nf-core pipeline using different datasets(ITS, PacBio, 18S, stingless bee microbiome,and IonTorrent data)
|        | Creating a test dataset, a test config and including flags for MBBU/16S-Accreditation QIIME2 pipeline
|        | Making a new documentation for MBBU/16S-Accreditation
| WEEK 6 | Running Yosef's pipeline
|        | Viewing nf-core/ampliseq functional analysis results using STAMP
| WEEK 7 | Report writing


## Pipelines' Functionality

* Some of the points in the criteria will be shown in reference to this range: 
  * 1 - Very good
  * 2 - Good
  * 3 - Fairly good
  * 4 - Bad
  * 5 - Very bad
* It is important to note that the MBBU-16S_Accreditation-QIIME2 and MBBU-16S_Accreditation-DADA2 are in the same github repository.

| Criteria | nf-core | H3ABionet-SOPs | H3aBionet-TADA | mbbu-16S_Accreditation-QIIME2 | mbbu-16S_Accreditation-DADA2 |
| -------- | ------- | -------------- | -------------- | ----------------------------- | ---------------------------- |
| Sequences | 16S, 18S, ITS |	16S | 16S, ITS | 16S | 16S |
| Tools and Databases | [Tools](https://github.com/nf-core/ampliseq#pipeline-summary),[Databases](https://nf-co.re/ampliseq/parameters#taxonomic-database)| [Tools](https://github.com/h3abionet/H3ABionet-SOPs/blob/master/pages/genomics_analysis/16s-rRNA/16s-rRNA.md#tools-referred-to-in-sop-tools), [Databases](https://github.com/h3abionet/H3ABionet-SOPs/blob/master/pages/genomics_analysis/16s-rRNA/16s-rRNA.md#databases-referred-to-in-sop-databases) | Not well defined | [Tools](https://github.com/mbbu/16S_Accreditation/blob/main/Qiime2_report.md#qiime-nexflow-pipeline), [Databases](https://github.com/mbbu/16S_Accreditation/blob/main/Qiime2_Nextflow/modules/chimera.nf) | [Tools](https://github.com/mbbu/16S_Accreditation/blob/main/Dada2_report.md#set-up), [Databases](https://github.com/mbbu/16S_Accreditation/blob/main/Dada2_Pipeline/dada2_pipeline.R) |
| QIIME2 | Yes | Yes | No | Yes | No | 
| DADA2 | No | No | Yes | No | Yes |
| Practice Dataset and Metadata | [Available](https://github.com/nf-core/ampliseq/blob/master/conf/test.config) | [Available](https://github.com/h3abionet/H3ABionet-SOPs/blob/master/pages/genomics_analysis/16s-rRNA/16s-rRNA.md#h3abionet-assessment-exercises) | N/A | N/A | N/A |
| Versions | [9 Releases] (https://github.com/nf-core/ampliseq/releases) | N/A (It is an SOP) | [1 Release](https://github.com/h3abionet/TADA/releases) | [1 Release](https://github.com/mbbu/16S_Accreditation/releases) | [1 Release](https://github.com/mbbu/16S_Accreditation/releases) |
| Command Arguments | [Available](https://nf-co.re/ampliseq/parameters#taxonomic-database) | N/A (It is an SOP) | [Available](https://github.com/h3abionet/TADA/blob/master/docs/usage.md#full-list-of-arguments-example) | N/A | N/A |
| Results | [Available](https://nf-co.re/ampliseq/2.1.1/output) | N/A (It is an SOP) | N/A | [Available](https://github.com/mbbu/16S_Accreditation/tree/main/Results) | N/A (While running this R-script, some results are saved: [script](https://github.com/mbbu/16S_Accreditation/blob/main/Dada2_Pipeline/dada2_pipeline.R) | 
| Documentation | 1 | For the guidelines provided (2), Navigating to the SOP(4) | 3 | 4 | 4|
| Contributions and Support | [Guidelines and Slack](https://github.com/nf-core/ampliseq#contributions-and-support), 17 contributors, 75 stars, 50 forks | N/A | 7 contributors, 10 stars, 11 forks | 9 contributors, 1 star, 2 forks | 9 contributors, 1 star, 2 forks |
| Running on cloud | Yes [Results from AWS Cloud](https://nf-co.re/ampliseq/results#ampliseq/results-80b3cb8b05d3b596bd0a52866e7febe40ea497db/) | N/A (It is an SOP) | Yes [AWS configs](https://github.com/h3abionet/TADA/tree/master/conf) | N/A | N/A |
| Languages | [Codes](https://github.com/nf-core/ampliseq/search?l=Groovy&type=code) | N/A (It is an SOP) | [Codes](https://github.com/h3abionet/TADA/search?l=nextflow) | [Codes](https://github.com/mbbu/16S_Accreditation/search?l=html) | [Codes](https://github.com/mbbu/16S_Accreditation/search?l=html) |
| Issues | [366](https://github.com/nf-core/ampliseq/search?l=Groovy&type=issues) | N/A (It is an SOP) | [32](https://github.com/h3abionet/TADA/search?l=nextflow&type=issues) | [10](https://github.com/mbbu/16S_Accreditation/search?l=html&type=issues) | [10](https://github.com/mbbu/16S_Accreditation/search?l=html&type=issues) |
| Last updated | October 2021 | February 2019 | September 2021 | April 2021 | April 2021 |

## More information on individual pipelines
### H3ABionet-SOPs/16s-rRNA-1-0.html
* There are questions on operation,run-time and output analysis that one would consider having as criteria in reviewing workflows.
* It would be a good SOP for an individual intending to create a QIIME pipeline from scratch.

## h3abionet/TADA
* This pipeline is a Targeted Amplicon Diversity Analysis(TADA) using DADA2, implemented in Nextflow.
* The typical command for running the pipeline is:
```
nextflow run uct-cbio/16S-rDNA-dada2-pipeline --reads '*_R{1,2}.fastq.gz' --trimFor 24 --trimRev 25 --reference 'gg_13_8_train_set_97.fa.gz' -profile uct_hex
```
* It outputs results mostly in RDS format.

## nf-core/ampliseq pipeline
* It supports paired-end Illumina or single-end Illumina, PacBio and IonTorrent data. 
* Analysis of 16S rRNA gene amplicons sequenced paired-end with Illumina is the default analysis.
* It runs with Conda, Docker, Podman,Shifter, Charliecloud or Singularity.
* Command for running the pipeline:
  ```
  nextflow run nf-core/ampliseq -profile <docker/singularity/podman/shifter/charliecloud/conda/institute> --input "path/to/data" --FW_primer "forward-primer-sequence" --RV_primer "reverse-primer-sequence" --metadata "Path/to/metadata/file""
  ```

![Image of how it runs and output expected](https://github.com/nf-core/ampliseq/blob/master/docs/images/ampliseq_workflow.png)


## MBBU 16S-Accreditation
### DADA2 MBBU 16S-Accreditation pipeline

* Steps in the pipeline:
![DADA2 MBBU 16S-Accreditation pipeline](https://github.com/mbbu/Reviewing-16s-Analysis-Workflows/blob/main/MBBU-16S-Accreditation-Dada2-Pipeline-Steps%20(2).png)

### QIIME MBBU 16S-Accreditation Pipeline

* It is summarized as follows:
![image](https://user-images.githubusercontent.com/91982522/149777440-1efe7a27-8034-492e-944d-d9edaa7b35ed.png)

## h3abionet16S
* We had difficulties running this and due to lack of an easy set-up, we did away with it.

## Workflow Comparison
This workflow comparison is only for the pipelines that we tested and found to be running without difficulties or with minimal difficulties.

<table>
    <tr>
       <th>Criteria</th>
       <th>nf-core/ampliseq</th>
       <th>TADA</th>
       <th>MBBU Accreditation Qiime</th>
       <th>MBBU Accreditation dada2</th>
   </tr>
   <tr>
       <td>Runtime</td>
       <td>21 minutes</td>
       <td>46 minutes</td>
       <td>#</td>
       <td>#</td>
   </tr>
   <tr>
       <td>Setup</td>
       <td>Easy, 1 command</td>
       <td>Easy, 1 command</td>
       <td>Hard, edit configs</td>
       <td>Hard</td>
   </tr>
   <tr>
       <td>Documentation</td>
       <td>Well documented</td>
       <td>Well documented</td>
       <td>Lacks setup instructions</td>
       <td>Not well documented</td>
   </tr>
    <tr>
       <td>Gaps</td>
       <td>None</td>
       <td>Test data, visualization</td>
       <td>Test config, functional analysis</td>
       <td>Not automated</td>
  
    
</table>

## Challenges
### Time
We had many holidays and a December break that in turn lessened the time to work on our objectives

### Internet
There was no internet for a period of time and could not progress with our work and had to change location of work

### Server space
Some processes would not run because of having less space in the server

## Recommendations
We purposed to do the following but hope some
* Nextflow automation of MBBU/16S-Accreditation DADA2 pipeline
* Adding picrust to MBBU/16S-Accreditation DADA2 pipeline for functional analysis
