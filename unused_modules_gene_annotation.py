# -*- coding: utf-8 -*-
"""
Created on Mon Oct  8 11:40:28 2018

@author:  sdeblank
"""

import requests
import re
import json
import vcf as pyvcf

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

#def ENSEMBLE_gene():
#stable_id="ENSG00000141510"
#request=requests.get('http://rest.ensembl.org/overlap/id/'+ stable_id +'?feature=transcript;content-type=application/json')
#response=request.text
#x=json.loads(response)
#print(x)


#with open("/home/cog/sdeblank/Downloads/homo_sapiens.json","r") as file:
#    content=file.readline()
#    json_parsed=json.loads(content)
#    print (json_parsed)
#    print (json_parsed["assembly_name"])
#    print (json_parsed["assembly_date"])
#
#    for i in range(0,len(json_parsed["top_level_region"])):
#        if (json_parsed["top_level_region"][i]["coord_system"] == "chromosome"):
#            print (json_parsed["top_level_region"][i]["name"] , "\t" , json_parsed["top_level_region"][i]["length"])
#
#server = "http://rest.ensembl.org"
#endpoint = "/xrefs/symbol/homo_sapiens/BRCA2"
#headers = {}
#headers['Content-Type'] = 'application/json'
#
#request=requests.get(server + endpoint, headers=headers)
#response=request.text
#data=json.loads(response)
#print (data)
#
#def find_gene_ID(gene):
#
#    server="http://rest.ensembl.org"
#    endpoint="/xrefs/symbol/homo_sapiens/"+gene
#    headers = {}
#    headers['Content-Type'] = 'application/json'
#
#    request=requests.get(server + endpoint, headers=headers)
#    response=request.text
#    data=json.loads(response)
#
#    return (data[0]["id"])
#
#def find_transcripts(ensemble_id):
#    server="http://rest.ensembl.org"
#    endpoint="/overlap/id/" + ensemble_id + "?feature=transcript"
#    headers = {}
#    headers['Content-Type'] = 'application/json'
#
#
#    request=requests.get(server + endpoint, headers=headers)
#    response=request.text
#    data=json.loads(response)
#
#    for entries in data:
#        return (str(entries['id'])+'\t'+"Location = " + str(entries['seq_region_name']) + ":" + str(entries['start']) + "-" + str(entries['end']) + "\t" + str(entries['biotype']) )
#
#    #print (data)
#
#
#searched_gene=input("Please type gene symbol: ")
#gene_ID=find_gene_ID(searched_gene)
#transcripts=find_transcripts(gene_ID)
