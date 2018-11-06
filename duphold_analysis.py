import vcf as pyvcf
import argparse

parser = argparse.ArgumentParser()
parser = argparse.ArgumentParser(description='Put here a description.')
parser.add_argument('vcf', help='VCF file')
parser.add_argument('-o', '--output', type=str, help='VCF output file', required=True)
args = parser.parse_args()

input=args.vcf
output=args.output

with open(input, "r") as input, open(output, "w") as output:
    VCF_READER=pyvcf.Reader(input)
    output.write("DHFC" + "\t" + "DHFFC" + "\t" + "DHBFC" + "\n")

    for record in VCF_READER:
        if isinstance(record.ALT[0], pyvcf.model._Breakend):
            if record.ALT[0].chr == record.CHROM and not record.ALT[0].orientation and record.ALT[0].remoteOrientation:
                print (record.ALT[0])
                DHFC=record.samples[0]["DHFC"]
                DHFFC=record.samples[0]["DHFFC"]
                DHBFC=record.samples[0]["DHBFC"]

                output.write(str(DHFC) + "\t" + str(DHFFC) + "\t" + str(DHBFC) + "\n")
