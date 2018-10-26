# -*- coding: utf-8 -*-
"""
Created on Mon Oct  8 14:59:36 2018

@author:  sdeblank
"""

import requests
import re
import json
import vcf as pyvcf
import sys
import argparse

parser = argparse.ArgumentParser()
parser = argparse.ArgumentParser(description='Put here a description.')
#parser.add_argument('vcf', help='VCF file')
parser.add_argument('-f', '--flank', default=200, type=int, help='Flank [default: 200]')
args = parser.parse_args()
FLANK=args.flank

def overlap_ENSEMBLE(bed_record):
    bed_record=bed_record.strip()
    bed_record=bed_record.split("\t")

    SERVER="https://rest.ensembl.org/overlap/region/"
    SPECIES="human"
    VERSION="GRCh37"
    CHROM=bed_record[0]
    SV_START=bed_record[1]
    SV_END=bed_record[2]
    ID=bed_record[3]
    HEADERS={"Content-Type" : "application/json"}


    region=bed_record[0] + ":" + bed_record[1] + "-" + bed_record[2]

    request = requests.get(SERVER+SPECIES+"/"+CHROM+":"+SV_START+"-"+SV_END+"?coord_system_version="+VERSION+";feature=gene", headers=HEADERS)
    response = request.text
    data=json.loads(response)
    genes={}
    if isinstance(data, list):
        genes[ID]=[gene["id"] for gene in data]
    else:
        print ("Error:", data["error"])
        genes[ID]="None"
    return genes

def alt_convert( record ):
    orientation = None
    remoteOrientation = None
    if record.INFO['SVTYPE'] == 'DEL':
        orientation = False
        remoteOrientation = True
    elif record.INFO['SVTYPE'] == 'DUP':
        orientation = True
        remoteOrientation = False
    elif 'INV3' in record.INFO:
    	orientation = False
    	remoteOrientation = False
    elif 'INV5' in record.INFO:
        orientation = True
        remoteOrientation = True

    if 'CT' in record.INFO:
        if record.INFO['CT'] == '3to5':
            orientation = False
            remoteOrientation = True
        elif record.INFO['CT'] == '3to3':
            orientation = False
            remoteOrientation = False
        elif record.INFO['CT'] == '5to3':
            orientation = True
            remoteOrientation = False
        elif record.INFO['CT'] == '5to5':
            orientation = True
            remoteOrientation = True
    if orientation is None or remoteOrientation is None:
        sys.exit()
    record.ALT = [ pyvcf.model._Breakend( record.CHROM, record.INFO['END'], orientation, remoteOrientation, record.REF, True ) ]
    return( record )


def bed_from_vcf(INPUT_VCF, OUTPUT_BED):
    with open(INPUT_VCF, "r") as vcf, open(OUTPUT_BED, "w") as OUTPUT_BED:
        VCF_READER=pyvcf.Reader(vcf)

        #OUTPUT_BED.write("CHROM" + "\t" + "Start_pos" + "\t" + "End_pos" + "\t" + "ID" + "\n")

        for record in VCF_READER:
            BEGIN_CHROM = str(record.CHROM)
            BEGIN_POS = record.POS
            ID=str(record.ID)



            if "INS" in record.ALT[0]:
                OUTPUT_BED.write(BEGIN_CHROM + "\t" + str(BEGIN_POS-FLANK) + "\t" + str(int(BEGIN_POS)+1+FLANK) + "\t"+ ID + "\n")

            elif not isinstance(record.ALT[0], pyvcf.model._Breakend):
                record = alt_convert(record)
            #If ALT is <INS> or something else --> Make ALT same as POS (Only location of insert needs to be taken into account, not the length)
            if  re.search("\w*[\[\]][\w\d]+:\d+[\[\]]\w*" ,str(record.ALT[0])) is not None:
                END_CHROM = record.ALT[0].chr #re.findall("[\w\d]+(?=:)", str(record.ALT[0]))[0]
                END_POS = record.ALT[0].pos #int(re.findall("(?<=:)[\w\d]+", str(record.ALT[0]))[0])

                #If ALT_chromosome != POS_chromosome --> Put both breakpoints on new line
                if END_CHROM != BEGIN_CHROM:
                    if record.ALT[0].orientation and record.ALT[0].remoteOrientation:
                        OUTPUT_BED.write(BEGIN_CHROM + "\t" + str(BEGIN_POS-1-FLANK) + "\t" + str(BEGIN_POS+FLANK) + "\t" + ID + "\n")
                        OUTPUT_BED.write(END_CHROM + "\t" + str(END_POS-1-FLANK) + "\t" + str(END_POS+FLANK) + "\t" + ID + "\n")
                        continue
                    elif not record.ALT[0].orientation and record.ALT[0].remoteOrientation:
                        OUTPUT_BED.write(BEGIN_CHROM + "\t" + str(BEGIN_POS-FLANK) + "\t" + str(BEGIN_POS+1+FLANK) + "\t" + ID + "\n")
                        OUTPUT_BED.write(END_CHROM + "\t" + str(END_POS-1-FLANK) + "\t" + str(END_POS+FLANK) + "\t" + ID + "\n")
                        continue
                    elif record.ALT[0].orientation and not record.ALT[0].remoteOrientation:
                        OUTPUT_BED.write(BEGIN_CHROM + "\t" + str(BEGIN_POS-1-FLANK) + "\t" + str(BEGIN_POS+FLANK) + "\t" + ID + "\n")
                        OUTPUT_BED.write(END_CHROM + "\t" + str(END_POS-FLANK) + "\t" + str(END_POS+1+FLANK) + "\t" + ID + "\n")
                        continue
                    elif not record.ALT[0].orientation and not record.ALT[0].remoteOrientation:
                        OUTPUT_BED.write(BEGIN_CHROM + "\t" + str(BEGIN_POS-FLANK) + "\t" + str(BEGIN_POS+1+FLANK) + "\t" + ID + "\n")
                        OUTPUT_BED.write(END_CHROM + "\t" + str(END_POS-FLANK) + "\t" + str(END_POS+1+FLANK) + "\t" + ID + "\n")
                        continue
                else:
                    OUTPUT_BED.write(BEGIN_CHROM + "\t" + str(BEGIN_POS-FLANK) + "\t" + str(END_POS+FLANK) + "\t" + ID + "\n")



            else:

                #sys.exit()
                continue

def annotate_vcf(INPUT_VCF, OUTPUT_VCF, ENSEMBLE_GENES):
    with open(INPUT_VCF, "r") as vcf, open(OUTPUT_VCF, "w") as OUTPUT_VCF:
        vcf.seek(0)
        VCF_READER=pyvcf.Reader(vcf)
        VCF_READER.infos['GENES']=pyvcf.parser._Info('GENES', ".", "String", "Overlapping genes", "NanoSV", "X")
        VCF_WRITER=pyvcf.Writer(OUTPUT_VCF, VCF_READER, lineterminator='\n')
        for record in VCF_READER:
            #record.INFO["GENES"]=",".join(ENSEMBLE_GENES[record.ID])
            if len(ENSEMBLE_GENES[record.ID])>0:
                record.add_info("GENES", ",".join(ENSEMBLE_GENES[record.ID]))
            VCF_WRITER.write_record(record)
            print (record.INFO)


#VCF_IN="/home/cog/sdeblank/Documents/sharc/COLO829_truth.vcf"
VCF_IN="/home/cog/sdeblank/Documents/sharc/test_short.vcf"
VCF_OUT=VCF_IN.replace(".vcf", "_gene_annotated.vcf")
BED=VCF_IN.replace(".vcf", ".bed")

bed_from_vcf(VCF_IN, BED)

with open(BED, "r") as BED_file:
    overlap_genes={}
    for entry in BED_file:
        overlap_genes.update(overlap_ENSEMBLE(entry))

annotate_vcf(VCF_IN, VCF_OUT, overlap_genes)
