

## Practical use case 1: inspecting a VCF file

We are going to follow the manual here https://eriqande.github.io/eca-bioinf-handbook/basic-handling-of-vcf-files.html

But before, we need to install bcftools. We will use ~~conda~~ mamba for this purpose (we have not explained what any of this is, we will, but for now we use mamba as it is equivalent and faster).

The first step will be to install mamba in the shell, for this purpose run the file `setup_micromamba.sh`.

```
bash setup_micromamba.sh
```

Now we will need to initialize micromamba:
```
microbamba activate
``` 

Now we can install bcftools
```
micromamba create --name bcftools bioconda::bcftools
```

And to activate the conda environment with bcfools:

```
micromamba activate bcftools
```

The data to use in the above manual is already downloaded and available in the `data` folder.


## Practical use case 2: Using a genomics viewer


### Step 0: install IGV

Install the Integrative Genomics Viewer (IGV) from:
https://software.broadinstitute.org/software/igv/download

Alternatively, you could use https://igv.org/app/ which does not require installation. But, they don't have exactly the same functinality and some of the questions will need some further work from your part.

In both cases you will need to download some of the files under https://github.com/priesgo/course-uas-intro-bioinformatics/tree/main/day_1/data

Open IGV.

### Step 1: load reference genome

Load the SARS-CoV-2 reference genome in FASTA format from `course-uas-intro-bioinformatics/day_1/practical_use_case_1/Sars_cov_2.ASM985889v3.dna.toplevel.fa` into IGV.

Optionally, download the reference from a public server.

`Genomes > Load genome from file... `

You should see now the single SARS-CoV-2 "chromosome" MN908947.3. The SARS-CoV-2 is almost 30,000 base pairs long. 

If you zoom in you will see the DNA letters that conform the reference genome. But where do genes start and end? This information is not in the FASTA file

### Step 2: load gene annotations

Load the gene annotations in GFF3 format from `course-uas-intro-bioinformatics/day_1/practical_use_case_1/gene_annotations.gff3` into IGV.

Optionally, download the gene annotations from a public server.

Right click on the new track and make sure that the visibility is `expanded`. You will now see all SARS-CoV-2 genes. The spike protein, the nucleocapsid protein, etc. If you zoom in you will now see in the new track the amino acids encoded in each gene.

### Step 3: load alignments

Load the test sample in BAM format from `course-uas-intro-bioinformatics/day_1/practical_use_case_1/alignments.bam` into IGV.

Now we should have two new tracks. On the top the coverage distribution showing the number of reads that overlap each position of the genome. On the bottom, the individual reads in the so called reads pile up.

You will observe that the reads have some positions colored. These are bases in the reads that differ from the reference genome. 

### Step 4: load the mutations file

Load the test sample in BAM format from `course-uas-intro-bioinformatics/day_1/practical_use_case_1/mutations.vcf` into IGV.

Now, we should have an additional track with the individual mutations mutations.


### Exercises

- What does the color coding of the mutations represent? (I suspect the color code may not be available in the web version...)
The color coding of the mutations represents the Variant Allele Frequency (VAF). A mutation with all reads supporting the mutation will have one unique color, a mutations with only half of the reads supporting it will have one half in one color and another half in another color, etc.

- Give me the VAF of mutation in position 3378?
- 0.0365854

- Which gene does the mutation at position 23403 affect? And which
- S

- What is the amino acid change and its position in the protein of mutation in genomic position 23403?
- D614G

- In position 17126 there is a deletion being called and it is supported only by two reads. Fetch the CIGAR string of the two reads supporting the deletion.
- 40M3D111M and 7S46M3D97M

- The genes in a virus are very close together, hence the stop codon of one gene and the start codon of the next are very close together. Take a screenshot with enough resolution that shows both stop and start codon.
