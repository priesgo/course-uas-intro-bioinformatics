# Practical use case 1: Using a genomics viewer

## Step 0: install prerequisites

If git is not installed already follow the instructions here https://git-scm.com/book/en/v2/Getting-Started-Installing-Git.

Clone the course repository:
```
git clone git@github.com:priesgo/course-uas-intro-bioinformatics.git
```
or 
```
git clone https://github.com/priesgo/course-uas-intro-bioinformatics.git
```

## Step 1: install IGV

Install the Integrative Genomics Viewer (IGV) from:
https://software.broadinstitute.org/software/igv/download

Open IGV.

## Step 2: load reference genome

Load the SARS-CoV-2 reference genome in FASTA format from `course-uas-intro-bioinformatics/day_1/practical_use_case_1/Sars_cov_2.ASM985889v3.dna.toplevel.fa` into IGV.

Optionally, download the reference from a public server.

`Genomes > Load genome from file... `

You should see now the single SARS-CoV-2 "chromosome" MN908947.3. The SARS-CoV-2 is almost 30,000 base pairs long. 

If you zoom in you will see the DNA letters that conform the reference genome. But where do genes start and end? This information is not in the FASTA file

## Step 3: load gene annotations

Load the gene annotations in GFF3 format from `course-uas-intro-bioinformatics/day_1/practical_use_case_1/gene_annotations.gff3` into IGV.

Optionally, download the gene annotations from a public server.

Right click on the new track and make sure that the visibility is `expanded`. You will now see all SARS-CoV-2 genes. The spike protein, the nucleocapsid protein, etc. If you zoom in you will now see in the new track the amino acids encoded in each gene.

## Step 4: load alignments

Load the test sample in BAM format from `course-uas-intro-bioinformatics/day_1/practical_use_case_1/alignments.bam` into IGV.

Now we should have two new tracks. On the top the coverage distribution showing the number of reads that overlap each position of the genome. On the bottom, the individual reads in the so called reads pile up.

You will observe that the reads have some positions colored. These are bases in the reads that differ from the reference genome. 

## Step 5: load the mutations file

Load the test sample in BAM format from `course-uas-intro-bioinformatics/day_1/practical_use_case_1/mutations.vcf` into IGV.

Now, we should have an additional track with the individual mutations mutations.


## Exercises

- What does the color coding of the mutations represent?
The color coding of the mutations represents the Variant Allele Frequency (VAF). A mutation with all reads supporting the mutation will have one unique color, a mutations with only half of the reads supporting it will have one half in one color and another half in another color, etc.

- Give me the VAF of mutation in position 3378?
- 0.0365854

- Which gene does the mutation 23403 affect?
- S

- What is the amino acid change and its position of mutation in genomic position 23403?
- D614G

- In position 17126 there is a deletion being called and it is supported only by two read. Fetch the CIGAR string of the two reads supporting the deletion.
- 40M3D111M and 7S46M3D97M