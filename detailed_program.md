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
**Practical use case 1**: manipulating a VCF file with bcftools
**Practical use case 2**: using a genomic brower, Integrative Genomics Viewer (IGV)


## Day 2

### Session 1
Version control with git
**Demo 1**: using git
Online biological databases and genome browsers
**Demo 2**: using biological databases
**Demo 3**: online genome browsers

### Session 2
Computational genomic pipelines - part II
Replicate and reproduce
Evidences of the reproducibility crisis
FAIR principles
Open science, open source

### Session 3 - hands on
**Practical use case 1**: git exercises following https://github.com/eficode-academy/git-katas
**Practical use case 2**: gathering data from a genomic database using a REST API


## Day 3

### Session 1
Software engineering best practices
Definition of legacy code
Automation and continuous integration
**Demo 1**: continuous integration environment
Software registries in bioinformatics (ie: PyPI, CRAN, bioconductor, bioconda)
**Demo 2**: using and managing conda environments

### Session 2
Programming in the biomedical sciences
Python crash course
Python libraries to handle popular bioinformatics data formats (pysam, cyvcf2, pybedtools)
Python libraries for machine learning (sklearn)
**Demo 3**: training a classification model
Artificial Intelligence (AI) tools for development. Copilot and others
**Demo 4**: assisted development by AI

### Session 3
**Practical use case 1**: train a machine learning model
**Practical use case 2**: develop a variant calling pipeline in bash
**Practical use case 3** (optional): develop a (simplistic) variant annotator in Python


## Day 4

### Session 1
Introduction to workflow managers or how to write better pipelines
Common features of a workflow manager
Python, R, Java or any of them?
Existing alternatives (eg: Common Workflow Language (CWL), Workflow Description Language (WDL), Galaxy, SnakeMake, Apache Airflow, Nextflow)

### Session 2
The Nextflow universe and the NF-core community
Crash course into Nextflow
Processes, subworkflows, channels and the DAG, (Directed Acyclic Graph)
Integration via configuration with different computational environments
Managing environments with conda and containerization with docker/singularity

### Session 3
**Practical use case 1**: migrate the bash pipeline into Nextflow
**Practical use case 2**: automate execution of the pipeline into the GitHub Continuous Integration environment