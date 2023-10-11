## Practical use case 1: develop a variant calling pipeline in bash

We are going to start from our test SARS-CoV-2 FASTQs and hopefully end with a usable set of mutations in a VCF file. 
For this purpose we need to use a number of bioinformatic tools that we will install with micromamba, we will process the data in the command line 
and finally to make it reusable we will write it down in a bash script.

### Step 1: alignment

First, create the environment:
```
micromamba create -f src/alignment.yaml
```

The file `alignment.yaml` contains all the dependencies we need to perform the alignment: BWA 2 and samtools.

Activate the environment:
```
micromamba activate alignment
```

Now run the following:
```
bwa-mem2 mem ../day_1/data/Sars_cov_2.ASM985889v3.dna.toplevel.fa ../day_1/data/reads_1.fastq.gz ../day_1/data/reads_2.fastq.gz | \
    samtools view -uS - | \
    samtools sort - > exercise.bam
```

Confirm that your new BAM file has some alignments into it:
```
samtools view exercise.bam | less
```


### Step 2: variant calling

First, create the environment:
```
micromamba create -f src/variant_calling.yaml
```

The file `alignment.yaml` contains all the dependencies we need to perform the alignment: BWA 2 and samtools.

Activate the environment:
```
micromamba activate variant_calling
```

Now run the following:
```
lofreq call  \
    --ref ../day_1/data/Sars_cov_2.ASM985889v3.dna.toplevel.fa \
    --call-indels \
    <( lofreq indelqual --dindel --ref ../day_1/data/Sars_cov_2.ASM985889v3.dna.toplevel.fa exercise.bam ) | bgzip > exercise.vcf.gz

bcftools index exercise.vcf.gz
```

Check that your VCF contains mutations:
```
bcftools view exercise.vcf.gz
```


### Step 3: putting it all together in a script

There is a template of a bash script in `src/script.sh`. 
We want to avoid writing in the script paths to files so it is reusable for different samples, for that reason it uses input parameters.

Fill the script with the alignment and variant calling commands using the defined parameters.


## Practical use case 2: develop a simple variant annotator in python

Our aim now is to implement some custom code in Python that reads a VCF file and adds an 
arbitrary annotation to every mutation.

We will use `cyvcf2` library to parse VCF files and `pybedtools` to parse BED files.
We will also make a quick overview of the `pandas` library used to manage data tables.

First, create the environment with these dependencies:
```
micromamba create -f src/annotation.yaml
```

The file `annotation.yaml` contains all the dependencies we need.

Activate the environment:
```
micromamba activate annotation
```

This will be the skeleton of our annotation that reads a VCF file, iterates through every mutation and writes it back to another VCF `src/annotation.py`.

```
from cyvcf2 import VCF, Writer, Variant

vcf = VCF('exercise.vcf.gz')
vcf_writer = Writer('annotated.vcf', vcf)

for variant in vcf:
    # it should add some info here to the variant
    vcf_writer.write_record(variant)

vcf_writer.close()
```

Now we want to read a BED file with some annotations and we put this in a convenient pandas dataframe.
This BED file contains the SARS-CoV-2 protein domains.

```
import pybedtools

bed = pybedtools.BedTool('data/pfam_names.bed.gz')
intersect = bed.intersect('exercise.vcf.gz')
df = intersect.to_dataframe()
```

Given the dataframe and a genomic coordinate we can perform a match. This will be our annotation.
```
match = df[(df.start <= variant.POS) & (df.end >= variant.POS)].name
if match.shape[0] > 0:
    variant.INFO["pfam"] = match.iloc[0]
```


If there is a match from the BED file we will put this into an INFO annotation.
Beware that the description of the new INFO annotation has to be added to the VCF beforehand.
```
from cyvcf2 import VCF, Writer, Variant
import pybedtools


bed = pybedtools.BedTool('data/pfam_names.bed.gz')
intersect = bed.intersect('exercise.vcf.gz')
df = intersect.to_dataframe()


vcf = VCF('exercise.vcf.gz')
vcf.add_info_to_header({'ID': 'pfam', 'Description':'', 'Type':'String', 'Number':1})
vcf_writer = Writer('annotated.vcf', vcf)


for variant in vcf:
    match = df[(df.start <= variant.POS) & (df.end >= variant.POS)].name
    if match.shape[0] > 0:
        variant.INFO["pfam"] = match.iloc[0]
    vcf_writer.write_record(variant)

vcf_writer.close()
```

Put it all together in the `src/annotation.py` and run it through one of our examples.
