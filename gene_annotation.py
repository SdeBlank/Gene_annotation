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

#############################################   CREATE A LIST OF CANCER GENES FROM THE TCGA DATABASE   #############################################

def create_TCGA_gene_list(CANCER_TYPE, MIN_SUPPORT):

############################ FILTER OUT ALL TCGA PROJECTS THAT WERE USED IN ANALYSING THE PRIMARY SITE CANCER (THUS FILTERING OUT PROJECTS LIKE FM-AD)
    SERVER_PROJECTS="https://api.gdc.cancer.gov/projects"
    FILTERS_PROJECTS={"op":"AND","content":[
                                {"op":"in","content":{"field":"projects.primary_site","value":CANCER_TYPE}},
                                {"op":"in","content":{"field":"projects.program.name","value":"TCGA"}}
                                ]}
    PARAMS_PROJECTS = {
        "filters": json.dumps(FILTERS_PROJECTS),
        "format": "JSON",
        "size": "1000"
        }

    request_projects=requests.get(SERVER_PROJECTS, params=PARAMS_PROJECTS)
    response_projects=request_projects.text
    response_projects=json.loads(response_projects)
    hits_projects=response_projects["data"]["hits"]
    PROJECTS=[]
    for project in hits_projects:
        PROJECTS.append(project["project_id"])
    print ("TCGA projects associated with "+CANCER_TYPE+": "+ ", ".join(PROJECTS))

############################# FILTER OUT ALL CASES THAT WERE USED IN THESE PROJECTS

    SERVER_CASES="https://api.gdc.cancer.gov/cases"
    FILTERS_CASES={"op":"AND","content":[
                            {"op":"in","content":{"field":"cases.primary_site","value":CANCER_TYPE}},
                            {"op":"in","content":{"field":"genes.is_cancer_gene_census","value":"true"}},
                            {"op":"in","content":{"field":"cases.project.project_id","value":PROJECTS}}
                            ]}
    PARAMS_CASES = {
        "filters": json.dumps(FILTERS_CASES),
        "format": "JSON",
        "size": "25000"
        }

    request_cases=requests.get(SERVER_CASES, params=PARAMS_CASES)
    response_cases=request_cases.text
    response_cases=json.loads(response_cases)
    hits_cases=response_cases["data"]["hits"]
    CASES=[]
    for hit in hits_cases:
        CASES.append(hit["submitter_id"])
    print ("Total number of cases", len(CASES))

############################ FILTER OUT ALL CASES THAT WERE USED IN MUTATION ANALYSIS
    SERVER_CASETYPE="https://api.gdc.cancer.gov/cases"

    CASE_NUMBER=0
    slice_start=0
    slice_end=0

    while slice_end < len(CASES):
        slice_end+=300
        if slice_end > len(CASES):
            slice_end=len(CASES)
        FILTERS_CASETYPE={"op":"in","content":{"field":"submitter_id","value":CASES[slice_start:slice_end]}}
        PARAMS_CASETYPE = {
            "filters": json.dumps(FILTERS_CASETYPE),
            "format": "JSON",
            "expand": "files",
            "size": "300"
            }
        request_casetype=requests.get(SERVER_CASETYPE, params=PARAMS_CASETYPE)
        response_casetype=request_casetype.text
        response_casetype=json.loads(response_casetype)
        hits_casetype=response_casetype["data"]["hits"]
        for hit in hits_casetype:
            for files in hit["files"]:
                file_type=files['data_category']
                if file_type == "Simple Nucleotide Variation":
                    CASE_NUMBER+=1
                    break

        slice_start+=300

    print ("Number of cases used in analysis", CASE_NUMBER)

############################ FILTER OUT ALL TCGA GENES THAT HAVE A GIVEN OCCURRENCE PERCENTAGE

    SERVER_GENES="https://api.gdc.cancer.gov/analysis/top_mutated_genes_by_project"
    FILTERS_GENES={"op":"AND","content":[
                                {"op":"in","content":{"field":"case.primary_site","value":CANCER_TYPE}}    ,
                                {"op":"in","content":{"field":"cases.project.project_id","value":PROJECTS}},
                                {"op":"in","content":{"field":"genes.is_cancer_gene_census","value":"true"}}
                                ]}
    FIELDS_GENES = [
        "gene_id",
        "symbol"
        ]

    FIELDS_GENES = ",".join(FIELDS_GENES)

    PARAMS_GENES = {
        "filters": json.dumps(FILTERS_GENES),
        "fields": FIELDS_GENES,
        "format": "JSON",
        "size": "100000"
        }

    request_genes=requests.get(SERVER_GENES, params=PARAMS_GENES)
    response_genes=request_genes.text
    response_genes=json.loads(response_genes)
    hits_genes=response_genes['data']["hits"]

    SIGNIFICANT_GENES={}
    for HIT in hits_genes:
        OCCURRENCE=float(float(HIT["_score"])/float(CASE_NUMBER))
        if OCCURRENCE>=float(MIN_SUPPORT):
            SIGNIFICANT_GENES[HIT["gene_id"]]=OCCURRENCE
    print (len(SIGNIFICANT_GENES))
    print ("Selecting genes with a minimal occurrence of "+str(MIN_SUPPORT)+"/"+str(CASE_NUMBER)+"="+str(float(MIN_SUPPORT)*CASE_NUMBER))
    return SIGNIFICANT_GENES

    # SIGNIFICANT_GENES=[]
    # for HIT in hits_genes:
    #     if float(HIT["_score"])/float(CASE_NUMBER)>=float(MIN_SUPPORT):
    #         SIGNIFICANT_GENES.append(HIT["gene_id"])
    #
    # print ("Selecting genes with a minimal occurence of "+str(MIN_SUPPORT)+"/"+str(CASE_NUMBER)+"="+str(float(MIN_SUPPORT)*CASE_NUMBER))
    # return SIGNIFICANT_GENES

#############################################   OVERLAP GENES THAT OVERLAP WITH GIVEN SV VCF AND TCGA CANCER GENES   #############################################
def vcf_annotate_tcga_genes_overlap(INPUT_VCF, OUTPUT_VCF, PROS_GENES, REGIONS):
    with open(INPUT_VCF, "r") as INPUT, open (OUTPUT_VCF, "w") as OUTPUT:
        count_pros_overlap=0
        count_gene_no_overlap=0
        count_no_gene=0
        x=0
        VCF_READER=pyvcf.Reader(INPUT)
        VCF_READER.infos['TCGAGENES']=pyvcf.parser._Info('TCGAGENES', 1, "Integer", "Number of prostate cancer genes overlapping with the SV region (+"+str(FLANK)+"bp flanking region)", "NanoSV", "X")
        VCF_WRITER=pyvcf.Writer(OUTPUT, VCF_READER, lineterminator='\n')
        GO={}
        for record in VCF_READER:
            x+=1
            overlap=0

            ENSEMBLE_OVERLAP=REGIONS[record.ID]["GENES"]
            OVERLAP=set(ENSEMBLE_OVERLAP).intersection(PROS_GENES)
            record.INFO["TCGAGENES"]=len(OVERLAP)
            VCF_WRITER.write_record(record)
            if len(OVERLAP)>0:
                score=0
                for i in OVERLAP:
                    #print (round(KNOWN_GENES[i], 3))
                    score+=(KNOWN_GENES[i] ** 2)
                score=int(round(score*100000, 0))
                if "SVLEN" in record.INFO:
                    print (str(record.ID) + "\t LENGTH=" + str(record.INFO["SVLEN"][0]) + "\t" + str(record.ALT[0]) + "\t" + "TCGAGENES=" + str(len(OVERLAP)) + "\t" + "SCORE=" + str(score))
                else:
                    print (str(record.ID) + "\t" + "TRANS/INS" + "\t" + "TCGAGENES=" + str(len(OVERLAP)) + "\t" + "SCORE=" + str(score))
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

#REGIONS=regions_from_vcf(VCF_IN)
#OVERLAP=overlap_ENSEMBLE(REGIONS)
KNOWN_GENES=create_TCGA_gene_list(CANCERTYPE, MIN_SUPPORT)
#vcf_annotate_tcga_genes_overlap(VCF_IN, VCF_GENE_SELECTED, KNOWN_GENES, OVERLAP)
