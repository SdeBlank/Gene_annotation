# -*- coding: utf-8 -*-
"""
Created on Mon Oct  8 14:59:36 2018

@author:  sdeblank
"""

import requests
import json
import vcf as pyvcf
import sys
import argparse

parser = argparse.ArgumentParser()
parser = argparse.ArgumentParser(description='Put here a description.')
#parser.add_argument('vcf', help='VCF file')
parser.add_argument('-f', '--flank', default=200, type=int, help='Flank [default: 200]')
parser.add_argument('-g', '--genes', default=None, help='Specific list of genes to compare VCF calls with [default: None]')

args = parser.parse_args()
FLANK=args.flank

def overlap_ENSEMBLE(bed_record):
    bed_record=bed_record.strip()
    bed_record=bed_record.split("\t")

    SERVER="https://GRCh37.rest.ensembl.org/overlap/region/"
    SPECIES="human"
    #VERSION="GRCh37"
    CHROM=bed_record[0]
    SV_START=int(bed_record[1])
    SV_END=int(bed_record[2])
    ID=bed_record[3]
    HEADERS={"Content-Type" : "application/json"}
    genes={}

    region=bed_record[0] + ":" + bed_record[1] + "-" + bed_record[2]

    if SV_END-SV_START <= 5000000:
        request = requests.get(SERVER+SPECIES+"/"+CHROM+":"+str(SV_START)+"-"+str(SV_END)+"?feature=gene", headers=HEADERS)
        response = request.text
        data=json.loads(response)
        if isinstance(data, list):
            genes[ID]=[gene["id"] for gene in data]
        else:
            print ("Error:", data["error"])
            genes[ID]=[]
        return genes
    else:
        TEMP_SV_START=SV_START
        TEMP_SV_END=TEMP_SV_START
        genes[ID]=[]

        while TEMP_SV_END < SV_END:
            TEMP_SV_END+=4999999
            if TEMP_SV_END > SV_END:
                TEMP_SV_END=SV_END
            request = requests.get(SERVER+SPECIES+"/"+CHROM+":"+str(TEMP_SV_START)+"-"+str(TEMP_SV_END)+"?feature=gene", headers=HEADERS)
            response = request.text
            data=json.loads(response)
            if isinstance(data, list):
                genes[ID]=genes[ID]+[gene["id"] for gene in data]
            else:
                print ("Error:", data["error"])
                genes[ID]=[]
                break
            TEMP_SV_START=TEMP_SV_END
        return genes

def gene_ontology(ENSEMBLE_ID):
    SERVER="https://GRCh37.rest.ensembl.org/xrefs/id/"
    HEADERS={"Content-Type" : "application/json"}

    request = requests.get(SERVER+ENSEMBLE_ID+"?external_db=GO;all_levels=1", headers=HEADERS)
    response = request.text
    data=json.loads(response)

    GO={}
    if isinstance(data, list):
        GO[ENSEMBLE_ID]=list(set([GO["description"] for GO in data]))
    else:
        print ("Error:", data["error"])
        GO[ENSEMBLE_ID]=[]

    return GO

def alt_convert( record ):
    orientation = None
    remoteOrientation = None
    if record.INFO['SVTYPE'] == 'DEL':
        orientation = False
        remoteOrientation = True
    elif record.INFO['SVTYPE'] == 'DUP':
        orientation = True
        remoteOrientation = False
    elif record.INFO['SVTYPE'] == 'TRA':
        strands=record.INFO['STRANDS']
        if strands == "++":
            orientation = False
            remoteOrientation = False
        elif strands == "+-":
            orientation = False
            remoteOrientation = True
        elif strands == "-+":
            orientation = True
            remoteOrientation = False
        elif strands == "--":
            orientation = True
            remoteOrientation = True
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

        for record in VCF_READER:
            BEGIN_CHROM = str(record.CHROM)
            BEGIN_POS = record.POS
            ID=str(record.ID)

            if "INS" in str(record.ALT[0]):
                OUTPUT_BED.write(BEGIN_CHROM + "\t" + str(BEGIN_POS-FLANK) + "\t" + str(int(BEGIN_POS)+1+FLANK) + "\t"+ ID + "\n")

            else:
                if not isinstance(record.ALT[0], pyvcf.model._Breakend):
                    record = alt_convert(record)

                END_CHROM = record.ALT[0].chr
                END_POS = record.ALT[0].pos

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


def annotate_vcf(INPUT_VCF, OUTPUT_VCF, ENSEMBLE_GENES):
    with open(INPUT_VCF, "r") as vcf, open(OUTPUT_VCF, "w") as OUTPUT_VCF:
        vcf.seek(0)
        VCF_READER=pyvcf.Reader(vcf)
        VCF_READER.infos['GENES']=pyvcf.parser._Info('GENES', ".", "String", "Genes overlapping with the SV region (+"+str(FLANK)+"bp flanking region)", "NanoSV", "X")
        VCF_WRITER=pyvcf.Writer(OUTPUT_VCF, VCF_READER, lineterminator='\n')
        for record in VCF_READER:
            #record.INFO["GENES"]=",".join(ENSEMBLE_GENES[record.ID])
            if len(ENSEMBLE_GENES[record.ID])>0:
                record.add_info("GENES", ",".join(ENSEMBLE_GENES[record.ID]))
            VCF_WRITER.write_record(record)

def create_gene_list(JSON):
    with open (JSON, "r") as JSON:
        PROS_GENES=[]
        for gene in JSON:
            PROS_GENES.append(gene['gene_id'])


#VCF_IN="/home/cog/sdeblank/Documents/sharc/SURVIVOR.vcf"
VCF_IN="/home/cog/sdeblank/Documents/sharc/EMC026T.nanosv.SHARC.primers.vcf"
#VCF_IN="/home/cog/sdeblank/Documents/sharc/test.vcf"
VCF_OUT=VCF_IN.replace(".vcf", "_gene_annotated.vcf")
BED=VCF_IN.replace(".vcf", ".bed")
KNOWN_GENES="/home/cog/sdeblank/Documents/sharc/pros_genes_min3.json"

bed_from_vcf(VCF_IN, BED)

with open(BED, "r") as BED_file:
    overlap_genes={}
    for entry in BED_file:
        overlap_genes.update(overlap_ENSEMBLE(entry))

annotate_vcf(VCF_IN, VCF_OUT, overlap_genes)
with open(VCF_OUT, "r") as vcf:
    count1=0
    count2=0
    count3=0
    VCF_READER=pyvcf.Reader(vcf)
    GO={}
    for record in VCF_READER:
        overlap=0
        if "GENES" in record.INFO:
            GENES=record.INFO["GENES"]
            for GENE in GENES:
                GO.update(gene_ontology(GENE))
            with open (KNOWN_GENES, "r") as pros:
                pros=json.loads(pros.read())
                PROS_GENES=[]
                for gene in pros:
                    PROS_GENES.append(gene['gene_id'])
                overlap=set(record.INFO['GENES']).intersection(PROS_GENES)
                if len(overlap)>0:
                    if "SVLEN" in record.INFO:
                        print (str(record.ID) + "\t LENGTH=" + str(record.INFO["SVLEN"][0]) + "\t" + str(record.ALT[0]) + "\t" + str(len(overlap)))
                    else:
                        print (str(record.ID) + "\t" + str(record.CHROM)+ "\t" + str(record.ALT[0]) + "\t" + str(len(overlap)))
                    count1+=1
                else:
                    count2+=1
        else:
            count3 +=1

    # count=0
    # for ENS_ID, GO_LIST in GO.items():
    #     for GO_TERM in GO_LIST:
    #         if GO_TERM=="nucleus":
    #              count += 1

    print ("Totaal: "+ str(x))
    print ("Gene overlap: ",count1)
    print ("Wel gen, geen overlap: ",count2)
    print ("Geen gen: ", count3)
