#!/usr/bin/python

import vcf as pyvcf
import sys
import argparse

parser = argparse.ArgumentParser()
parser = argparse.ArgumentParser(description='Put here a description.')
parser.add_argument('-v', '--vcf', type=str, help='VCF file', required=True)
parser.add_argument('-o', '--output', type=str, help='Output', required=True)
args = parser.parse_args()

def select_sv_length(vcf):
    print ("ID" + "\t" + "SV length")
    with open(vcf, "r") as vcf:
        VCF_READER=pyvcf.Reader(vcf)

        for record in VCF_READER:
            if "SVLEN" in record.INFO:
                print(str(record.ID) + "\t" + str(record.INFO["SVLEN"][0]))
            else:
                print(str(record.ID) + "\t" + "0")


SHARC_primer_vcf=args.vcf
select_sv_length(SHARC_primer_vcf)
