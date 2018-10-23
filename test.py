import requests
import re
import json
import vcf as pyvcf

#print(float(int(7)/int(3)))
CANCER_TYPE=["TCGA-OV", "TCGA-SARC"]

# SERVER_CASES="https://api.gdc.cancer.gov/cases"
SERVER_CASES="https://api.gdc.cancer.gov/analysis/mutated_cases_count_by_project"
SERVER_GENES="https://api.gdc.cancer.gov/analysis/top_mutated_genes_by_project"

FIELDS_GENES = [
    "gene_id",
    "_score",
    "symbol"
    ]

FIELDS_GENES = ",".join(FIELDS_GENES)

FIELDS_CASES = [
    "project_id"
    ]

FIELDS_CASES = ",".join(FIELDS_CASES)

#FILTERS_CASES={"op":"AND","content":[{"op":"in","content":{"field":"cases.primary_site","value":["Prostate gland"]}}]}
# FILTERS_CASES={"op":"in","content":{"field":"primary_site","value":"Prostate gland"}}
FILTERS_CASES={"op":"AND","content":[{"op":"in","content":{"field":"project.project_id","value":"TCGA-OV"}}]}
FILTERS_GENES={"op":"AND","content":[{"op":"in","content":{"field":"case.project.project_id","value":[CANCER_TYPE]}}]}

PARAMS_GENES = {
    "filters": json.dumps(FILTERS_GENES),
    "fields": FIELDS_GENES,
    "format": "JSON",
    "size": "10"
    }

PARAMS_CASES = {
    "filters": json.dumps(FILTERS_CASES),
    "format": "JSON",
    "size": "100"
    }

request_cases=requests.get(SERVER_CASES, params=PARAMS_CASES)
response_cases=request_cases.text
response_cases=json.loads(response_cases)
hits_cases=response_cases["aggregations"]["projects"]["buckets"]
CASE_NUMBER=0
for case in hits_cases:
    CASE_NUMBER+=int(hits_cases[0]["case_summary"]["case_with_ssm"]["doc_count"])
#print(hits_cases)
#print (len(hits_cases))
# CASE_NUMBER=int(hits_cases[0]['summary']['case_count'])
print(CASE_NUMBER)

request_genes=requests.get(SERVER_GENES, params=PARAMS_GENES)
response_genes=request_genes.text
response_genes=json.loads(response_genes)
#print (response_genes)
hits_genes=response_genes['data']["hits"]


# a=0
# SIGNIFICANT_GENES=[]
# for HIT in hits:
#     print (HIT)
#     if int(HIT["_score"])>=5:
#         a+=1
#         SIGNIFICANT_GENES.append(HIT)
# print (a)


# VCF_IN="/home/cog/sdeblank/Documents/sharc/EMC026T.nanosv.SHARC.primers.vcf"
#
# with open(VCF_IN, "r") as vcf, open("/home/cog/sdeblank/Documents/sharc/SV_length.vcf", "w") as OUTPUT_VCF:
#     vcf.seek(0)
#     VCF_READER=pyvcf.Reader(vcf)
#     VCF_WRITER=pyvcf.Writer(OUTPUT_VCF, VCF_READER, lineterminator='\n')
#     for record in VCF_READER:
#         #record.INFO["GENES"]=",".join(ENSEMBLE_GENES[record.ID])
#         if "SVLEN" in record.INFO:
#             if record.INFO["SVLEN"][0] > 1000000:
#                 VCF_WRITER.write_record(record)

# dic={'A':1, 'B':1000, 'C':5, 'D':4}
# sorted_dic=[(k, dic[k]) for k in sorted(dic, key=dic.get, reverse=True)]
# print (sorted_dic)
# print (dic.items())
# for key, value in sorted_dic:
#     print (key + ": " + str(value))
#
# dic={}
# x = []
# y = [3,5,8]
# dic["a"]=[]
# dic["a"]=dic["a"]+y
# z=x+y
#
# print (dic["a"])
#
# for x in range(5):
#     print ("x", x)
#     for y in range(10):
#         print ("y", y)
#         if y==6:
#             break
# print ("Done")
