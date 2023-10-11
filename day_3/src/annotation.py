#!/env/python
from cyvcf2 import VCF, Writer, Variant
import argparse
import logging


parser = argparse.ArgumentParser(description="variant annotator",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     epilog="Copyright (c) 2023 UAS")
parser.add_argument("--input-vcf", dest="input_vcf", action="store", help="The VCF to annotate", required=True)
parser.add_argument("--output-vcf", dest="output_vcf", action="store", help="The annotated VCF", required=True)


args = parser.parse_args()

logging.info("annotator starting")



# read input VCF
# vcf = VCF('exercise.vcf.gz')
vcf = VCF(args.input_vcf)

# create a writer to a new VCF
# vcf_writer = Writer('annotated.vcf', vcf)
vcf_writer = Writer(args.output_vcf, vcf)


# iterates through every variant
for variant in vcf:

    # TODO: it should add some info here to the variant
    
    # write the annotated to the output VCF
    vcf_writer.write_record(variant)


# close the output VCF
vcf_writer.close()