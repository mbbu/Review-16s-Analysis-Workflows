# Targeted Amplicon Diversity Analysis Workflow 
Assignees: Nelly and Collins

### Background
The 16S ribosomal RNA is necessary for the synthesis of all prokaryotic proteins.

Reproducibility is very vital in research as it ensures that results are the same no matter how many times a computational pipeline is executed. Workflow languages such as Nexflow and containers such as Docker or Singularity are key tools in ensuring reproducibility in research.

### Aim
This mini-project aims to review existing microbiome workflows, identify great ones and extend the workflows where there are gaps, especially to make them useful in insect and pathogen data. Review the workflows using the following criteria:
1.  How easy are they to set up and use? Do they provide accessible documentation and tutorials?
2. Are they fast and easily scalable based on available compute resources?
3. Can they scale to the cloud? 
4. Can they be used on a variety of data, including insects and pathogen microbiome
5. Are they implemented in the latest specifications and versions of the tools? For example, Whether the pipeline implements Nextflow DSL2 syntax and docker or singularity containers
6. Are they well and regularly maintained? When were they updated last?

### Some existing pipelines
- https://github.com/nf-core/ampliseq 
- https://h3abionet.github.io/H3ABionet-SOPs/16s-rRNA-1-0.html
- https://github.com/h3abionet/TADA
- https://github.com/mbbu/16S_Accreditation
- https://github.com/nanoporetech/pipeline-transcriptome-de/tree/paired_dge_dtu
-  https://github.com/biocorecrg/master_of_pores
-  https://github.com/h3abionet/h3abionet16S

### Tasks 
1. Create a Roadmap for the miniproject.
2. Test existing 16srRNA pipelines. 
3. Address the identified gap. 
4. Document your work clearly on GitHub using wikis and GitHub pages.
5. Document the papers you are reading, a link to the paper, and a sentence or two on why you included them.


**N/B:**
- Use test dataset from either of the listed pipelines.
- Demonstrate collaborative research skills, informative visualization, and report writing

## Some useful resources and references 
1. [16S rRNA, One of the Most Important rRNAs](https://www.cd-genomics.com/blog/16s-rrna-one-of-the-most-important-rrnas/)
2. [QIIME Tutorials¶](http://qiime.org/tutorials/)
3. [DADA2 Pipeline Tutorial (1.8)](https://benjjneb.github.io/dada2/tutorial_1_8.html)
4. [Microbiome Helper: a Custom and Streamlined Workflow for Microbiome Research](https://journals.asm.org/doi/10.1128/mSystems.00127-16)
5. 