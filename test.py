import requests
import re
import json
import vcf as pyvcf

SERVER="https://api.gdc.cancer.gov/projects"

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
                "field":"project_id",
                "value":[
                    "TCGA-PRAD"
                ]
            }
        }
    ]
}

PARAMS = {
    "filters": json.dumps(FILTERS),
#    "fields": FIELDS,
    "format": "JSON",
    "size": "5"
    }

request=requests.get(SERVER, params=PARAMS)
response=request.text
response=json.loads(response)
hits=response['data']["hits"]
print (response)
#hits=json['data']["hits"]

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
