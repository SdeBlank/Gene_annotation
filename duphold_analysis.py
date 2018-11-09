import vcf as pyvcf
import argparse

parser = argparse.ArgumentParser()
parser = argparse.ArgumentParser(description='Put here a description.')
parser.add_argument('vcf', help='VCF file')
parser.add_argument('-o', '--output', type=str, help='VCF output file', required=True)
args = parser.parse_args()

input=args.vcf
output=args.output

def distribution(input, output):
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

def check_real(input, output):
    with open(input, "r") as input, open(output, "w") as output:
        VCF_READER=pyvcf.Reader(input)

        for record in VCF_READER:
            if isinstance(record.ALT[0], pyvcf.model._Breakend):
                if record.ALT[0].chr == record.CHROM:
                    DHFC=str(record.samples[0]["DHFC"])
                    DHFFC=str(record.samples[0]["DHFFC"])
                    DHBFC=str(record.samples[0]["DHBFC"])

                    if DHFFC=="inf":
                        DHFFC="10"
                    elif "nan" in DHFFC:
                        DHFFC="0"
                    if not record.ALT[0].orientation and record.ALT[0].remoteOrientation and float(DHFFC)<=0.5:
                        output.write(str(record.ID)+ "\t" + "DEL" + "\t" + DHFC+":"+DHFFC+":"+DHBFC + "\t" + str(record.FILTER) + "\n")
                    elif record.ALT[0].orientation and not record.ALT[0].remoteOrientation and float(DHFFC)>=1.5:
                        output.write(str(record.ID)+ "\t" + "DUP" + "\t" + DHFC+":"+DHFFC+":"+DHBFC + "\t" + str(record.FILTER) + "\n")
                    elif (record.ALT[0].orientation and record.ALT[0].remoteOrientation) or ( not record.ALT[0].orientation and not record.ALT[0].remoteOrientation):
                        output.write(str(record.ID)+ "\t" + "INV" + "\t" + DHFC+":"+DHFFC+":"+DHBFC + "\t" + str(record.FILTER) + "\n")
                    else:
                        print(str(record.ID))
                else:
                    output.write(str(record.ID)+ "\t" + "TRA" + "\t" + "NA" + "\t" + str(record.FILTER) + "\n")
            else:
                output.write(str(record.ID)+ "\t" + str(record.ALT[0]) + "\t" + "NA" + "\t" + str(record.FILTER) + "\n")

check_real(input, output)
