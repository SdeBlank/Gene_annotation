import requests
import json
import vcf as pyvcf
import sys
import argparse

def rank_GO(gene_JSON):
    with open(gene_JSON, "r") as JSON:
        JSON=JSON.read()
        JSON=json.loads(JSON)

    GO_count={}
    for entry in JSON:
        gene_id=entry['gene_id']
        SERVER="https://rest.ensembl.org/xrefs/id/"
        #SERVER="https://GRCh37.rest.ensembl.org/xrefs/id/"
        HEADERS={"Content-Type" : "application/json"}

        request = requests.get(SERVER+gene_id+"?external_db=GO;all_levels=1", headers=HEADERS)
        response = request.text
        data=json.loads(response)

        GO={}
        if isinstance(data, list):
            GO[gene_id]=list(set([(GO["display_id"], GO["description"]) for GO in data]))
        else:
            print ("Error:", data["error"])
            #genes[ID]="None"


        for id in GO:
            GO_LIST=GO[id]
            #print (GO_LIST)
            for description in GO_LIST:
                if description not in GO_count:
                    GO_count[description]=1
                else:
                    GO_count[description]+=1

    ranked_GO=[(ID, GO_count[ID]) for ID in sorted(GO_count, key=GO_count.get, reverse=True)]
    for id, count in ranked_GO:
        if count >=10:
            print (str(id) + ": " + str(count))


file="/home/cog/sdeblank/Documents/sharc/genes.json"
rank_GO(file)
