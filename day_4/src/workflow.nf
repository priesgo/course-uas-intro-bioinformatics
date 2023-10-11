#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


// these are the parameters to the workflow
params.fastq1 = false
params.fastq2 = false
params.fasta = false
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
    # TODO: fill in the call to the bwa command
    """
}

process VARIANT_CALLING {
    cpus 1
    memory '2g'
    publishDir "${output}", mode: "copy", pattern: "*.bcf"
    
    input:
        val(fasta)
        path(bam)
        val(output)

    output:
        file("${output}.bcf")

    """
    # TODO: fill in the call to the lofreq command
    """
}


workflow {

    ALIGNMENT(params.fasta, file(params.fastq1), file(params.fastq2), params.output)
    VARIANT_CALLING(params.fasta, ALIGNMENT.out, params.output)

    //PFAM_ANNOTATOR(VARIANT_CALLING.out.vcf)

}