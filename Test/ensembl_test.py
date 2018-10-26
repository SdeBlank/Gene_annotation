# -*- coding: utf-8 -*-
"""
Created on Mon Oct  8 11:40:28 2018

@author:  sdeblank
"""

import requests
import re
import json
import vcf as pyvcf



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

