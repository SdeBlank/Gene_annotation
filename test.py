import requests
import re
import json
import vcf as pyvcf
import sys

file1="/data/sharc/truth_ids"
file2="/data/sharc/colo829.full.ranked"
file3="/data/sharc/colo829.full.ranked.truth"

with open(file1, "r") as file1, open(file2, "r") as file2, open(file3, "w") as output:
    truth=[]
    for id in file1:
        id=id.strip()
        truth.append(str(id))
    for line in file2:
        columns=line.strip().split("\t")
        if str(columns[0]) in truth:
            print ("JA")
            columns.append("TRUTH")
            output.write(str("\t".join(columns))+"\n")
        else:
            output.write(line)



# TRY=1
# while TRY < 10:
#     print (TRY)
#     if TRY==5:
#         Y=5
#     try:
#         X=Y
#         print(X)
#         break
#     except:
#         print ("NO")
#         if TRY==7:
#             sys.exit("Error while requesting from ICGC database after 3 tries")
#         TRY +=1
#
# print ("BLALBLABLA")
# x=0
# for n in li:
#     if n>x:
#         x=n
#         a.insert(0, n)
#     elif n>=3:
#         a.append(n)
# print (a)

# SERVER_CASES="https://api.gdc.cancer.gov/cases"
# FILTERS_CASES={"op":"AND","content":[
#                         {"op":"in","content":{"field":"primary_site","value":"Ovary"}},
#                         {"op":"in","content":{"field":"project.project_id","value":["TCGA-OV", "TCGA-SARC"]}}
#                         ]}
# PARAMS_CASES = {
#     "filters": json.dumps(FILTERS_CASES),
#     "format": "JSON",
#     "size": "25000"
#     }
#
# request_cases=requests.get(SERVER_CASES, params=PARAMS_CASES)
# response_cases=request_cases.text
# response_cases=json.loads(response_cases)
# hits_cases=response_cases["data"]["hits"]
# CASES=[]
# for hit in hits_cases:
#     CASES.append(hit["submitter_id"])
# print ("Total number of cases", len(CASES))
#
# SERVER_CASES="https://api.gdc.cancer.gov/ssm_occurrences"
#
# CASE_NUMBER=0
# for case in CASES:
#     FILTERS_CASES={"op":"AND","content":[
#                                     {"op":"in","content":{"field":"case.submitter_id","value":case}}    #,
#                                     #{"op":"in","content":{"field":"cases.project.project_id","value":["TCGA-OV", "TCGA-SARC"]}}
#                                     ]}
#
#     PARAMS_CASES = {
#         "filters": json.dumps(FILTERS_CASES),
#         "format": "JSON",
#         "size": "1"
#         }
#
#     request_cases=requests.get(SERVER_CASES, params=PARAMS_CASES)
#     response_cases=request_cases.text
#     response_cases=json.loads(response_cases)
#     hits_cases=response_cases["data"]["hits"]
#     print (hits_cases)
#     if len(hits_cases)>0:
#         CASE_NUMBER+=1
#
# print (CASE_NUMBER)



#     hits_cases=response_cases["aggregations"]["projects"]["buckets"]
#     CASE_NUMBER=0
#     for case in hits_cases:
#         if
#         CASE_NUMBER+=int(hits_cases[0]["case_summary"]["case_with_ssm"]["doc_count"])
# #print(hits_cases)
# #print (len(hits_cases))
# # CASE_NUMBER=int(hits_cases[0]['summary']['case_count'])
# print(CASE_NUMBER)


#
# ###################################################################################
# SERVER_CASETYPE="https://api.gdc.cancer.gov/cases"
#
# CASE_NUMBER=0
# slice_start=0
# slice_end=0
#
# while slice_end < len(CASES):
#     slice_end+=300
#     if slice_end > len(CASES):
#         slice_end=len(CASES)
#     FILTERS_CASETYPE={"op":"in","content":{"field":"submitter_id","value":CASES[slice_start:slice_end]}}
#     PARAMS_CASETYPE = {
#         "filters": json.dumps(FILTERS_CASETYPE),
#         "format": "JSON",
#         "expand": "files",
#         "size": "300"
#         }
#     request_casetype=requests.get(SERVER_CASETYPE, params=PARAMS_CASETYPE)
#     response_casetype=request_casetype.text
#     response_casetype=json.loads(response_casetype)
#     hits_casetype=response_casetype["data"]["hits"]
#     for hit in hits_casetype:
#         for files in hit["files"]:
#             file_type=files['data_category']
#             print (file_type)
#             if file_type == "Simple Nucleotide Variation":
#                 CASE_NUMBER+=1
#                 break
#
#     slice_start+=300
#     print (slice_end)
#
# print ("Number of cases with SNVs", CASE_NUMBER)
# CASE_NUMBER=0
# for count, case in enumerate(CASES):             ####################################   #REMOVE COUNT!!!
#     print (count+1)
#     FILTERS_CASETYPE={"op":"in","content":{"field":"submitter_id","value":case}}
#     PARAMS_CASETYPE = {
#         "filters": json.dumps(FILTERS_CASETYPE),
#         "format": "JSON",
#         "expand": "files",
#         "size": "100"
#         }
#     request_casetype=requests.get(SERVER_CASETYPE, params=PARAMS_CASETYPE)
#     response_casetype=request_casetype.text
#     response_casetype=json.loads(response_casetype)
#     hits_casetype=response_casetype["data"]["hits"][0]["files"]
#     for files in hits_casetype:
#         file_type=files['data_category']
#         if file_type == "Simple Nucleotide Variation":
#             CASE_NUMBER+=1
#             break
# print ("Number of cases with SNVs", CASE_NUMBER)

# SERVER_CASES="https://api.gdc.cancer.gov/cases"
# FILTERS_CASES={"op":"AND","content":[
#                         {"op":"in","content":{"field":"primary_site","value":"Ovary"}},
#                         {"op":"in","content":{"field":"project.project_id","value":["TCGA-OV", "TCGA-SARC"]}}
#                         ]}
#                         #{"op":"in","content":{"field":"project.project_id","value":["TCGA-OV", "TCGA_SARC"]}}
# PARAMS_CASES = {
#     "filters": json.dumps(FILTERS_CASES),
#     "format": "JSON",
#     "size": "25000"
#     }
#
# request_cases=requests.get(SERVER_CASES, params=PARAMS_CASES)
# response_cases=request_cases.text
# response_cases=json.loads(response_cases)
# hits_cases=response_cases["data"]["hits"]
# CASES=[]
# for hit in hits_cases:
#     CASES.append(hit["submitter_id"])
# print ("Total number of cases", len(CASES))

# #print(float(int(7)/int(3)))
# CANCER_TYPE=["TCGA-OV", "TCGA-SARC"]
#

#
# request_genes=requests.get(SERVER_GENES, params=PARAMS_GENES)
# response_genes=request_genes.text
# response_genes=json.loads(response_genes)
# print (response_genes)
# hits_genes=response_genes['data']["hits"]
# #print (hits_genes)
#
# a=0
# SIGNIFICANT_GENES=[]
# for HIT in hits_genes:
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
