## Practical use case 1: develop a Nextflow variant calling workflow

We are going now to implement a similar pipeline to what we implemented on day 3, but this time in Nextflow.

Nextflow is a workflow manager that allows to integrate different technologies seamlessly in the same workflow 
and that allows to implement good practices in terms of reproducibility.

### Step 0: install Nextflow

Go to your home folder `/home/your_user_name` and run `wget -qO- https://get.nextflow.io | bash`

Now move the nextflow binary to your bin folder so this is executable from anywhere.
```
mv nextflow bin/
```

### Step 1: write the alignment and variant calling processes

Open the workflow template `src/workflow.nf`. The parameters and the connection between processes of the 
workflow is implemented but you need to add in each separate process the BWA and lofreq command we were working with in bash.

Run the workflow as follows:
```
nextflow src/workflow.nf --fasta `pwd`/../day_1/data/Sars_cov_2.ASM985889v3.dna.toplevel.fa --fastq1 ../day_1/data/reads_1.fastq.gz --fastq2 ../day_1/data/reads_2.fastq.gz --output test
```

**NOTE**: beware that the fasta reference is passed as an absolute path. This is necessary so Nextflow finds the indexes 
that are sitting next to the FASTA file. Otherwise, Nextflow copies the input files to intermediate temporary folders.

You will need to activate the script environment we created yesterday: `micromamba activate script`

### Step 2: make use of isolated conda environments for each process

To reach a better reproducibility each process may have its own environment. This is especially important as the workflow 
and its dependencies grow.

Add the directive conda to each process with the needed depdendencies. For instance, for the alignment:
```
conda ("bioconda::bwa-mem2=2.2.1 bioconda::samtools=1.12")
``` 

Then you will need to edit the file `.nextflow.config` and set `conda.enabled = true`.

THIS IS FAILING WITH micromamba! => just keep using the large environment, not ideal... but works.


### Step 3: add a step to the workflow using the Python annotator

Create a subfolder `bin` under the folder `src`. Then copy your python script into that folder. 

You will be able to call the script from a Nextflow process.

Beware, that you will need to pass the pfam BED file to the workflow as a parameter.

**NOTE**: because step 2 is not working you will need to create a environment with the alignment and 
variant calling dependencies plus the dependencies from our custom annotator.

Put this into a Nextflow process and add this step to the workflow.


### Step 4: add a step that adds functional annotations to the variants

There are many implementation but we will use the one in `bcftools csq`, but beware that we need the latest version 
of bcftools 1.18. Update out bcftools environment to this version.

Then you will need to index your VCF file obtained yesterday: `tabix your.vcf` or you could use the VCF under `day_4/data/mutations.vcf.gz`.

Annotate with a command as follows:
```
bcftools csq --force -f ../day_1/data/Sars_cov_2.ASM985889v3.dna.toplevel.fa -g ../day_1/data/gene_annotations.gff3 data/mutations.vcf.gz -Ob -o data/mutations.csq.vcf.gz
``` 

Put this into a Nextflow process and add this step to the workflow.


### Step 5: test your workflow over a newly downloaded file

Go to https://www.covid19dataportal.org/

This URL in particular should point you to 2,346,868 SARS-CoV-2 samples with certain technical criteria that make them fit for our pipeline:
https://www.covid19dataportal.org/search/sequences?db=sra-experiment-covid19&size=15&facets=instrument_platform:ILLUMINA,library_strategy:AMPLICON,TAXON:2697049,host:Homo%20sapiens&crossReferencesOption=all

Click one randomly.

In the table under "Generated FASTQ files: FTP" select the files ending in "_1.fastq.gz" and "_2.fastq.gz" and click "get download script".

It should look like this:
```
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR211/033/SRR21180633/SRR21180633_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR211/033/SRR21180633/SRR21180633_2.fastq.gz
```

Run it in the cloud shell. And now run your workflow over it!

You can commit the resulting files to your repository.