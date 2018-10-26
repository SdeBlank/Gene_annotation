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
parser.add_argument('-s', '--support', default=0.01, type=float, help='Minimal percentage of cancer patients supporting the mutated gene [default: 0.01]')
parser.add_argument('-c', '--cancertype', type=str, help='Primary site of cancer', required=True)

args = parser.parse_args()

#############################################   CONVERT DIFFERENT VCF SV NOTATIONS TO bracket notations N[Chr:pos[   #############################################
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


#############################################   FILTER OUT REGIONS AROUND BREAKPOINTS FROM A VCF FILE   #############################################
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

#############################################   OVERLAP SV WITH KNOWN ENSEMBLE GENES   #############################################
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

#############################################   CREATE A LIST OF CANCER GENES FROM THE ICGC DATABASE   #############################################

def create_TCGA_gene_list(CANCER_TYPE, MIN_SUPPORT):

##################################### Top 1500 genes with occurence > given occurence
    SIGNIFICANT_GENES={}
    GENE_NAMES={}

    SERVER_GENES="https://dcc.icgc.org/api/v1/genes"
    FILTERS_GENES={
                    "donor":{"primarySite":{"is":[CANCER_TYPE]}
                    ,"availableDataTypes":{"is":["ssm"]}}
                    #,"mutation":{"functionalImpact":{"is":["High"]}}
                    #,"gene":{"curatedSetId":{"is":["GS1"]}}

                    }
    FILTERS_GENES=json.dumps(FILTERS_GENES)

    MAX_GENES=int(requests.get("https://dcc.icgc.org/api/v1/genes/count?filters="+FILTERS_GENES).text)
    CASE_NUMBER=int(requests.get("https://dcc.icgc.org/api/v1/donors/count?filters="+FILTERS_GENES).text)

    SLICE=0
    STOP=False
    print ("All",MAX_GENES, CASE_NUMBER)
    while SLICE < MAX_GENES and STOP==False :
        print(SLICE)
        PARAMS_GENES = {
            "filters": FILTERS_GENES,
            "format": "JSON",
            "size": "100",
            "from": str(SLICE)
            }
        request_genes=requests.get(SERVER_GENES, params=PARAMS_GENES)
        response_genes=request_genes.text
        response_genes=json.loads(response_genes)
        hits_genes=response_genes["hits"]

        for HIT in hits_genes:
            OCCURRENCE=float(float(HIT["affectedDonorCountFiltered"])/float(CASE_NUMBER))
            if OCCURRENCE>=float(MIN_SUPPORT):
                SIGNIFICANT_GENES[HIT["id"]]={}
                SIGNIFICANT_GENES[HIT["id"]]["score"]=1
                SIGNIFICANT_GENES[HIT["id"]]["symbol"]=HIT["symbol"]
            else:
                STOP=True
                break

        SLICE+=100

##################################### Genes with occurence > given occurence and cancer_census is true
    SERVER_GENES="https://dcc.icgc.org/api/v1/genes"
    FILTERS_GENES={
                    "donor":{"primarySite":{"is":[CANCER_TYPE]}
                    ,"availableDataTypes":{"is":["ssm"]}}
                    #,"mutation":{"functionalImpact":{"is":["High"]}}
                    ,"gene":{"curatedSetId":{"is":["GS1"]}}
                    }
    FILTERS_GENES=json.dumps(FILTERS_GENES)

    MAX_GENES=int(requests.get("https://dcc.icgc.org/api/v1/genes/count?filters="+FILTERS_GENES).text)
    CASE_NUMBER=int(requests.get("https://dcc.icgc.org/api/v1/donors/count?filters="+FILTERS_GENES).text)

    print ("Cancer_census",MAX_GENES, CASE_NUMBER)
    SLICE=0
    STOP=False
    while SLICE < MAX_GENES and STOP==False:
        print(SLICE)
        PARAMS_GENES = {
            "filters": FILTERS_GENES,
            "format": "JSON",
            "size": "100",
            "from": str(SLICE)
            }

        request_genes=requests.get(SERVER_GENES, params=PARAMS_GENES)
        response_genes=request_genes.text
        response_genes=json.loads(response_genes)
        hits_genes=response_genes["hits"]

        for HIT in hits_genes:
            OCCURRENCE=float(float(HIT["affectedDonorCountFiltered"])/float(CASE_NUMBER))
            if OCCURRENCE>=float(MIN_SUPPORT):
                SIGNIFICANT_GENES[HIT["id"]]["score"]=3
            else:
                STOP=True
                break
        SLICE+=100

##################################### Genes with occurence > 0.1 and cancer_census is true
    SERVER_GENES="https://dcc.icgc.org/api/v1/genes"
    FILTERS_GENES={
                    "donor":{"primarySite":{"is":[CANCER_TYPE]}
                    ,"availableDataTypes":{"is":["ssm"]}}
                    #,"mutation":{"functionalImpact":{"is":["High"]}}
                    ,"gene":{"curatedSetId":{"is":["GS1"]}}
                    }
    FILTERS_GENES=json.dumps(FILTERS_GENES)

    MAX_GENES=int(requests.get("https://dcc.icgc.org/api/v1/genes/count?filters="+FILTERS_GENES).text)
    CASE_NUMBER=int(requests.get("https://dcc.icgc.org/api/v1/donors/count?filters="+FILTERS_GENES).text)

    print ("Cancer_census + 0.1 occurrence",MAX_GENES, CASE_NUMBER)
    SLICE=0
    STOP=False
    while SLICE < MAX_GENES and STOP==False:
        print(SLICE)
        PARAMS_GENES = {
            "filters": FILTERS_GENES,
            "format": "JSON",
            "size": "100",
            "from": str(SLICE)
            }

        request_genes=requests.get(SERVER_GENES, params=PARAMS_GENES)
        response_genes=request_genes.text
        response_genes=json.loads(response_genes)
        hits_genes=response_genes["hits"]

        for HIT in hits_genes:
            OCCURRENCE=float(float(HIT["affectedDonorCountFiltered"])/float(CASE_NUMBER))
            if OCCURRENCE>=0.1:
                SIGNIFICANT_GENES[HIT["id"]]["score"]=4
            else:
                STOP=True
                break
        SLICE+=100

##################################### Genes that have known high impact mutations
    SERVER_GENES="https://dcc.icgc.org/api/v1/genes"
    FILTERS_GENES={
                    "donor":{"primarySite":{"is":[CANCER_TYPE]}
                    ,"availableDataTypes":{"is":["ssm"]}}
                    ,"mutation":{"functionalImpact":{"is":["High"]}}
                    #,"gene":{"curatedSetId":{"is":["GS1"]}}
                    }
    FILTERS_GENES=json.dumps(FILTERS_GENES)

    MAX_GENES=int(requests.get("https://dcc.icgc.org/api/v1/genes/count?filters="+FILTERS_GENES).text)
    CASE_NUMBER=int(requests.get("https://dcc.icgc.org/api/v1/donors/count?filters="+FILTERS_GENES).text)

    print ("Impact",MAX_GENES, CASE_NUMBER)
    SLICE=0
    while SLICE < MAX_GENES:
        PARAMS_GENES = {
            "filters": FILTERS_GENES,
            "format": "JSON",
            "size": "100",
            "from": str(SLICE)
            }

        request_genes=requests.get(SERVER_GENES, params=PARAMS_GENES)
        response_genes=request_genes.text
        response_genes=json.loads(response_genes)
        hits_genes=response_genes["hits"]

        for HIT in hits_genes:
            #OCCURRENCE=float(float(HIT["affectedDonorCountFiltered"])/float(CASE_NUMBER))
            if HIT["id"] in SIGNIFICANT_GENES:
                SIGNIFICANT_GENES[HIT["id"]]["score"]+=1
        SLICE+=100

##################################### Return
    print (len(SIGNIFICANT_GENES))
    print ("Selecting genes with a minimal occurrence of "+str(MIN_SUPPORT)+"/"+str(CASE_NUMBER)+"="+str(float(MIN_SUPPORT)*CASE_NUMBER))
    return SIGNIFICANT_GENES


#############################################   OVERLAP GENES THAT OVERLAP WITH GIVEN SV VCF AND TCGA CANCER GENES   #############################################
def vcf_annotate_tcga_genes_overlap(INPUT_VCF, OUTPUT_VCF, ICGC_GENES, REGIONS):
    with open(INPUT_VCF, "r") as INPUT, open (OUTPUT_VCF, "w") as OUTPUT:

        count_pros_overlap=0
        count_gene_no_overlap=0
        count_no_gene=0
        x=0
        VCF_READER=pyvcf.Reader(INPUT)
        VCF_READER.infos['ICGC_SCORE']=pyvcf.parser._Info('ICGC_SCORE', 1, "Integer", "Score of ICGC cancer genes overlapping with the SV region (+"+str(FLANK)+"bp flanking region)", "NanoSV", "X")
        VCF_READER.infos['ICGC_OVERLAP']=pyvcf.parser._Info('ICGC_SCORE', ".", "String", "ICGC cancer gene ordered from high ICGC_SCORE to low ICGC_SCORE overlapping with the SV region (+"+str(FLANK)+"bp flanking region)", "NanoSV", "X")
        # 0 = No overlap
        # 1 = Overlap ----> occurrene > given occurrence
        # 2 = Overlap ----> occurrene > given occurrence + impact==HIGH
        # 3 = Overlap ----> occurrene > given occurrence + cancer gene
        # 4 = Overlap ----> occurrene > given occurrence + impact==HIGH + cancer gene
        # 5 = Overlap ----> occurrene > 10% + impact==HIGH + cancer gene

        VCF_WRITER=pyvcf.Writer(OUTPUT, VCF_READER, lineterminator='\n')
        GO={}
        for record in VCF_READER:
            x+=1
            overlap=0

            ENSEMBLE_OVERLAP=REGIONS[record.ID]["GENES"]
            OVERLAP=set(ENSEMBLE_OVERLAP).intersection(ICGC_GENES)

            SCORE=0
            GENE=[]
            for gene in OVERLAP:
                gene_score=ICGC_GENES[gene]["score"]
                if gene_score>SCORE:
                    SCORE=gene_score
                    GENE.insert(0, ICGC_GENES[gene]["symbol"])
                elif gene_score>=3:
                    GENE.append(ICGC_GENES[gene]["symbol"])
            record.INFO["ICGC_SCORE"]=SCORE
            record.INFO["ICGC_OVERLAP"]=",".join(GENE)
            VCF_WRITER.write_record(record)
            if len(OVERLAP)>0:
                if "SVLEN" in record.INFO:
                    print (str(record.ID) + "\t LENGTH=" + str(record.INFO["SVLEN"][0]) + "\t" + str(record.ALT[0]) + "\t" + "ICGC_OVERLAP=" + str(len(OVERLAP)) + "\t" + "SCORE=" + str(SCORE))
                else:
                    print (str(record.ID) + "\t" + "TRANS/INS" + "\t" + "ICGC_OVERLAP=" + str(len(OVERLAP)) + "\t" + "SCORE=" + str(SCORE))
                count_pros_overlap+=1
            elif len(REGIONS[record.ID]["GENES"])>0:
                count_gene_no_overlap+=1
            elif len(REGIONS[record.ID]["GENES"])==0:
                count_no_gene +=1

        print ("\n" +"Totaal: "+ str(x))
        print ("Gene overlap: ",count_pros_overlap)
        print ("Wel gen, geen overlap: ",count_gene_no_overlap)
        print ("Geen gen: ", count_no_gene)


#############################################   RUNNING CODE   #############################################
FLANK=args.flank
MIN_SUPPORT=args.support
CANCERTYPE=args.cancertype

VCF_IN=args.vcf
VCF_GENE_SELECTED=VCF_IN.replace(".vcf", "_gene_selection.vcf")

REGIONS=regions_from_vcf(VCF_IN)
OVERLAP=overlap_ENSEMBLE(REGIONS)
KNOWN_GENES=create_TCGA_gene_list(CANCERTYPE, MIN_SUPPORT)
vcf_annotate_tcga_genes_overlap(VCF_IN, VCF_GENE_SELECTED, KNOWN_GENES, OVERLAP)
