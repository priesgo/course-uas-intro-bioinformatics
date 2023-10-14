#!/usr/bin/env python

from cyvcf2 import VCF, Writer, Variant
import argparse
import logging
import pybedtools


"""
TODO:
You will need to complete this script with the appropriate code to read
the BED file, match every mutation with potential overlapping protein domains 
from the BED file and finally add the annotation in the INFO field
"""


parser = argparse.ArgumentParser(description="variant annotator",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     epilog="Copyright (c) 2023 UAS")
parser.add_argument("--input-vcf", dest="input_vcf", action="store", help="The VCF to annotate", required=True)
parser.add_argument("--output-vcf", dest="output_vcf", action="store", help="The annotated VCF", required=True)
parser.add_argument("--bed", dest="bed", action="store", help="The BED file to annotate", required=True)


args = parser.parse_args()

logging.info("annotator starting")



# read input VCF
# vcf = VCF('exercise.vcf.gz')
vcf = VCF(args.input_vcf)

# we add to the header the new annotation that will be named "pfam"
vcf.add_info_to_header({'ID': 'pfam', 'Description':'', 'Type':'String', 'Number':1})

# create a writer to a new VCF
# vcf_writer = Writer('annotated.vcf', vcf)
vcf_writer = Writer(args.output_vcf, vcf)


bed = pybedtools.BedTool(args.bed)
intersect = bed.intersect(args.input_vcf)
df = intersect.to_dataframe()


# iterates through every variant
for variant in vcf:

    # TODO: it should add some annotation here to the variant
    match = df[(df.start <= variant.POS) & (df.end >= variant.POS)].name
    if match.shape[0] > 0:
        variant.INFO["pfam"] = match.iloc[0]

    # write the annotated to the output VCF
    vcf_writer.write_record(variant)


# close the output VCF
vcf_writer.close()