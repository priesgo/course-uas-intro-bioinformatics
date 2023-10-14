# Detailed program

## Day 1

### Session 1
Presentations
Definition of bioinformatics
Skills needed by a bioinformatician
Brief history of bioinformatics

### Session 2
Overview of NGS methods
Standard genomic data formats
Computational genomic pipelines - part I
**Demo 1**: manipulating a VCF file with bcftools
**Demo 2**: using a genomic brower, Integrative Genomics Viewer (IGV)

### Session 3 - hands on
Setup Google Cloud Shell and GitHub
Install micromamba
**Practical use case 2**: using a genomic brower, Integrative Genomics Viewer (IGV)


## Day 2

### Session 1
Version control with git
**Demo 1**: using git
**Practical use case 1**: git exercises following https://github.com/eficode-academy/git-katas
### Session 2
Online biological databases and genome browsers
**Demo 2**: using biological databases
**Demo 3**: online genome browsers
### Session 3 - hands on
**Practical use case 2**: manipulating a VCF file with bcftools
**Practical use case 3**: gathering data from a genomic database using a REST API


## Day 3
### Session 1

Replicate and reproduce
Evidences of the reproducibility crisis
FAIR principles
Open science, open source
### Session 2
Software engineering best practices
Definition of legacy code
Software registries in bioinformatics (ie: PyPI, CRAN, bioconductor, bioconda)
**Demo 2**: using and managing conda environments
Programming in the biomedical sciences
**Practical use case 1**: develop a variant calling pipeline in bash
### Session 3
Python crash course
Python libraries to handle popular bioinformatics data formats (pysam, cyvcf2, pybedtools)
**Practical use case 2**: develop a (simplistic) variant annotator in Python


## Day 4

### Session 1
Integration via configuration with different computational environments
Managing environments with conda and containerization with docker/singularity
Introduction to workflow managers or how to write better pipelines
Common features of a workflow manager
Python, R, Java or any of them?
Existing alternatives (eg: Common Workflow Language (CWL), Workflow Description Language (WDL), Galaxy, SnakeMake, Apache Airflow, Nextflow)
The Nextflow universe and the NF-core community
Crash course into Nextflow
Processes, subworkflows, channels and the DAG, (Directed Acyclic Graph)
**Practical use case 1**: migrate the bash pipeline into Nextflow

### Session 2
Computational genomic pipelines - part II
**Practical use case 2**: download public SARS-CoV-2 data from the COVID-19 data portal
### Session 3
**Practical use case 3**: add additional custom and functional annotations into the pipeline