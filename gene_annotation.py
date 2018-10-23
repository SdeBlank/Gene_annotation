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
parser.add_argument('vcf', help='VCF file')
parser.add_argument('-f', '--flank', default=200, type=int, help='Flank [default: 200]')
parser.add_argument('-s', '--support', default=200, type=int, help='Minimal supporting cases a gene has been found in a cancer type [default: 5]')

args = parser.parse_args()


def overlap_ENSEMBLE(REGIONS):
    SERVER="https://GRCh37.rest.ensembl.org/overlap/region/"
    SPECIES="human"
    HEADERS={"Content-Type" : "application/json"}

    for ID in REGIONS:
        REGIONS[ID]["GENES"]=[]

        for region in REGIONS[ID]["REGION"]:
            CHROM=region["Chrom"]
            SV_START=int(region["Start"])
            SV_END=int(region["End"])

            if SV_END-SV_START <= 5000000:
                request = requests.get(SERVER+SPECIES+"/"+CHROM+":"+str(SV_START)+"-"+str(SV_END)+"?feature=gene", headers=HEADERS)
                response = request.text
                data=json.loads(response)
                if isinstance(data, list):
                    REGIONS[ID]["GENES"]=REGIONS[ID]["GENES"]+[gene["id"] for gene in data]
                else:
                    print ("Error:", data["error"])
                    REGIONS[ID]["GENES"]=[]
            else:
                TEMP_SV_START=SV_START
                TEMP_SV_END=TEMP_SV_START
                while TEMP_SV_END < SV_END:
                    TEMP_SV_END+=4999999
                    if TEMP_SV_END > SV_END:
                        TEMP_SV_END=SV_END
                    request = requests.get(SERVER+SPECIES+"/"+CHROM+":"+str(TEMP_SV_START)+"-"+str(TEMP_SV_END)+"?feature=gene", headers=HEADERS)
                    response = request.text
                    data=json.loads(response)
                    if isinstance(data, list):
                        REGIONS[ID]["GENES"]=REGIONS[ID]["GENES"]+[gene["id"] for gene in data]
                    else:
                        print ("Error:", data["error"])
                        REGIONS[ID]["GENES"]=[]
                        break
                    TEMP_SV_START=TEMP_SV_END
    return REGIONS

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


def regions_from_vcf(INPUT_VCF):
    with open(INPUT_VCF, "r") as vcf:
        VCF_READER=pyvcf.Reader(vcf)

        SV_DATA={}

        for record in VCF_READER:
            BEGIN_CHROM = str(record.CHROM)
            BEGIN_POS = record.POS
            ID=str(record.ID)
            SV_DATA[ID]={"REGION":[]}

            if "INS" in str(record.ALT[0]):
                REGION_START=BEGIN_POS-FLANK
                REGION_END=BEGIN_POS+1+FLANK
                SV_DATA[ID]["REGION"].append({"Chrom":BEGIN_CHROM, "Start":REGION_START, "End":REGION_END})

            else:
                if not isinstance(record.ALT[0], pyvcf.model._Breakend):
                    record = alt_convert(record)

                END_CHROM = record.ALT[0].chr
                END_POS = record.ALT[0].pos

                if END_CHROM != BEGIN_CHROM:
                    if record.ALT[0].orientation and record.ALT[0].remoteOrientation:
                        REGION_START_1 = BEGIN_POS-1-FLANK
                        REGION_END_1 = BEGIN_POS+FLANK
                        SV_DATA[ID]["REGION"].append({"Chrom":BEGIN_CHROM, "Start":REGION_START_1, "End":REGION_END_1})
                        REGION_START_2 = END_POS-1-FLANK
                        REGION_END_2 = END_POS+FLANK
                        SV_DATA[ID]["REGION"].append({"Chrom":END_CHROM, "Start":REGION_START_2, "End":REGION_END_2})
                        continue
                    elif not record.ALT[0].orientation and record.ALT[0].remoteOrientation:
                        REGION_START_1 = BEGIN_POS-FLANK
                        REGION_END_1 = BEGIN_POS+1+FLANK
                        SV_DATA[ID]["REGION"].append({"Chrom":BEGIN_CHROM, "Start":REGION_START_1, "End":REGION_END_1})
                        REGION_START_2 = END_POS-1-FLANK
                        REGION_END_2 = END_POS+FLANK
                        SV_DATA[ID]["REGION"].append({"Chrom":END_CHROM, "Start":REGION_START_2, "End":REGION_END_2})
                        continue
                    elif record.ALT[0].orientation and not record.ALT[0].remoteOrientation:
                        REGION_START_1 = BEGIN_POS-1-FLANK
                        REGION_END_1 = BEGIN_POS+FLANK
                        SV_DATA[ID]["REGION"].append({"Chrom":BEGIN_CHROM, "Start":REGION_START_1, "End":REGION_END_1})
                        REGION_START_2 = END_POS-FLANK
                        REGION_END_2 = END_POS+1+FLANK
                        SV_DATA[ID]["REGION"].append({"Chrom":END_CHROM, "Start":REGION_START_2, "End":REGION_END_2})
                        continue
                    elif not record.ALT[0].orientation and not record.ALT[0].remoteOrientation:
                        REGION_START_1 = BEGIN_POS-FLANK
                        REGION_END_1 = BEGIN_POS+1+FLANK
                        SV_DATA[ID]["REGION"].append({"Chrom":BEGIN_CHROM, "Start":REGION_START_1, "End":REGION_END_1})
                        REGION_START_2 = END_POS-FLANK
                        REGION_END_2 = END_POS+1+FLANK
                        SV_DATA[ID]["REGION"].append({"Chrom":END_CHROM, "Start":REGION_START_2, "End":REGION_END_2})
                        continue
                else:
                    REGION_START = BEGIN_POS-FLANK
                    REGION_END = END_POS+FLANK
                    SV_DATA[ID]["REGION"].append({"Chrom":BEGIN_CHROM, "Start":REGION_START, "End":REGION_END})
        return (SV_DATA)

def annotate_genes_vcf(INPUT_VCF, OUTPUT_VCF, ENSEMBLE_GENES):
    with open(INPUT_VCF, "r") as vcf, open(OUTPUT_VCF, "w") as OUTPUT_VCF:
        vcf.seek(0)
        VCF_READER=pyvcf.Reader(vcf)
        VCF_READER.infos['GENES']=pyvcf.parser._Info('GENES', ".", "String", "Genes overlapping with the SV region (+"+str(FLANK)+"bp flanking region)", "NanoSV", "X")
        VCF_WRITER=pyvcf.Writer(OUTPUT_VCF, VCF_READER, lineterminator='\n')
        for record in VCF_READER:
            if len(ENSEMBLE_GENES[record.ID]["GENES"])>0:
                record.add_info("GENES", ",".join(ENSEMBLE_GENES[record.ID]["GENES"]))
            VCF_WRITER.write_record(record)

def create_gene_list(CANCER_TYPE, MIN_SUPPORT):

    SERVER="https://api.gdc.cancer.gov/analysis/top_mutated_genes_by_project"

    FIELDS = [
        "gene_id",
        "symbol"
        ]

    FIELDS = ",".join(FIELDS)

    FILTERS={
        "op":"AND",
        "content":[
            {
                "op":"in",
                "content":{
                    "field":"case.project.project_id",
                    "value":[
                        CANCER_TYPE
                    ]
                }
            }
        ]
    }

    PARAMS = {
        "filters": json.dumps(FILTERS),
        "fields": FIELDS,
        "format": "JSON",
        "size": "10000"
        }

    request=requests.get(SERVER, params=PARAMS)
    response=request.text
    response=json.loads(response)
    hits=response['data']["hits"]

    SIGNIFICANT_GENES=[]
    for HIT in hits:
        if int(HIT["_score"])>=MIN_SUPPORT:
            SIGNIFICANT_GENES.append(HIT["gene_id"])
    return SIGNIFICANT_GENES

def vcf_annotate_pros_genes_overlap(INPUT_VCF, OUTPUT_VCF, PROS_GENES, REGIONS):
    with open(INPUT_VCF, "r") as INPUT, open (OUTPUT_VCF, "w") as OUTPUT:
        count_pros_overlap=0
        count_gene_no_overlap=0
        count_no_gene=0
        x=0
        VCF_READER=pyvcf.Reader(INPUT)
        VCF_READER.infos['PROSGENE']=pyvcf.parser._Info('PROSGENE', 1, "Integer", "Number of prostate cancer genes overlapping with the SV region (+"+str(FLANK)+"bp flanking region)", "NanoSV", "X")
        VCF_WRITER=pyvcf.Writer(OUTPUT, VCF_READER, lineterminator='\n')
        GO={}
        for record in VCF_READER:
            x+=1
            overlap=0

            ENSEMBLE_OVERLAP=REGIONS[record.ID]["GENES"]
            OVERLAP=set(ENSEMBLE_OVERLAP).intersection(PROS_GENES)
            record.INFO["PROSGENE"]=len(OVERLAP)
            VCF_WRITER.write_record(record)
            if len(OVERLAP)>0:
                #record.add_info("PROSGENE", len(OVERLAP))
                if "SVLEN" in record.INFO:
                    print (str(record.ID) + "\t LENGTH=" + str(record.INFO["SVLEN"][0]) + "\t" + str(record.ALT[0]) + "\t" + "PROSGENES=" + str(len(OVERLAP)))
                else:
                    print (str(record.ID) + "\t" + "TRANS/INS" + "\t" + "PROSGENES=" + str(len(OVERLAP)))
                count_pros_overlap+=1
            elif len(REGIONS[record.ID]["GENES"])>0:
                count_gene_no_overlap+=1
            elif len(REGIONS[record.ID]["GENES"])==0:
                count_no_gene +=1

        print ("\n" +"Totaal: "+ str(x))
        print ("Gene overlap: ",count_pros_overlap)
        print ("Wel gen, geen overlap: ",count_gene_no_overlap)
        print ("Geen gen: ", count_no_gene)

FLANK=args.flank
MIN_SUPPORT=int(args.support)

#KNOWN_GENES="/home/cog/sdeblank/Documents/sharc/pros_genes_min5.json"
VCF_IN=args.vcf
#VCF_IN="/home/cog/sdeblank/Documents/sharc/SURVIVOR.vcf"
#VCF_IN="/home/cog/sdeblank/Documents/sharc/EMC026T.nanosv.SHARC.primers.vcf"


#VCF_ANNOTATED=VCF_IN.replace(".vcf", "_gene_annotated.vcf")
VCF_GENE_SELECTED=VCF_IN.replace(".vcf", "_gene_selection.vcf")
#BED=VCF_IN.replace(".vcf", ".bed")

REGIONS=regions_from_vcf(VCF_IN)
KNOWN_GENES=create_gene_list("TCGA-PRAD", MIN_SUPPORT)
OVERLAP=overlap_ENSEMBLE(REGIONS)
vcf_annotate_pros_genes_overlap(VCF_IN, VCF_GENE_SELECTED, KNOWN_GENES, OVERLAP)
