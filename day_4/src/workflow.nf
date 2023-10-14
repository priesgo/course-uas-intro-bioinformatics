#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


// these are the parameters to the workflow
params.fastq1 = false
params.fastq2 = false
params.fasta = false
params.bed = false
params.gff = false
params.output = false

if (! params.fasta) {
    log.error "--fasta is required"
    exit 1
}
if (! params.fastq1) {
    log.error "--fastq1 is required"
    exit 1
}
if (! params.fastq2) {
    log.error "--fastq2 is required"
    exit 1
}
if (! params.output) {
    log.error "--output is required"
    exit 1
}
if (! params.bed) {
    log.error "--bed is required"
    exit 1
}
if (! params.gff) {
    log.error "--gff is required"
    exit 1
}

/*
this is a comment
*/

// this is a comment
process ALIGNMENT {
    cpus 1
    memory '2g'
    publishDir "${output}", mode: "copy", pattern: "*.bam"

    input:
        val(fasta)
        path(fastq1)
        path(fastq2)
        val(output)

    output:
        file("${output}.bam")

    """
    bwa-mem2 mem $fasta \
        $fastq1 $fastq2 | \
        samtools view -uS - | \
        samtools sort - > ${output}.bam
    """
}

process VARIANT_CALLING {
    cpus 1
    memory '2g'
    publishDir "${output}", mode: "copy", pattern: "*.vcf.gz"
    
    input:
        val(fasta)
        path(bam)
        val(output)

    output:
        file("${output}.vcf.gz")

    """
    lofreq call  \
        --ref $fasta \
        --call-indels \
        <( lofreq indelqual --dindel --ref $fasta $bam ) | bgzip > ${output}.vcf.gz

    bcftools index ${output}.vcf.gz
    """
}


process PFAM_ANNOTATOR {
    cpus 1
    memory '2g'
    publishDir "${output}", mode: "copy", pattern: "*.vcf.gz"
    
    input:
        path(bed)
        path(vcf)
        val(output)

    output:
        file("${output}.pfam.vcf.gz")

    """
    annotation.py --input-vcf ${vcf} \
        --output-vcf ${output}.pfam.vcf.gz \
        --bed ${bed}
    """
}

process FUNCTIONAL_ANNOTATOR {
    cpus 1
    memory '2g'
    publishDir "${output}", mode: "copy", pattern: "*.vcf.gz"
    
    input:
        val(fasta)
        val(gff)
        path(vcf)
        val(output)

    output:
        file("${output}.csq.vcf.gz")

    """
    bcftools csq --force \
        -f ${fasta} \
        -g ${gff} \
        ${vcf} \
        -Ob -o ${output}.csq.vcf.gz
    """
}


workflow {

    ALIGNMENT(params.fasta, file(params.fastq1), file(params.fastq2), params.output)
    VARIANT_CALLING(params.fasta, ALIGNMENT.out, params.output)    
    PFAM_ANNOTATOR(file(params.bed), VARIANT_CALLING.out, params.output)
    FUNCTIONAL_ANNOTATOR(params.fasta, params.gff, PFAM_ANNOTATOR.out, params.output)

}