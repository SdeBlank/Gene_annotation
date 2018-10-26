# -*- coding: utf-8 -*-
"""
Created on Mon Oct  8 14:59:36 2018

@author:  sdeblank
"""

import requests
import re
import json
import vcf as pyvcf


# with open("/home/cog/sdeblank/Documents/sharc/test.vcf", "r") as file:
#     vcf=pyvcf.Reader(file)
#     print(vcf.infos)
#     # print ("CHROM" + "\t" + "Start_pos" + "\t" + "End_pos" + "\t" + "ID")
#     for record in vcf:
#         if not isinstance(record.ALT[0], pyvcf.model._Breakend):
#             print (record.ALT[0])

    #     BEGIN_CHROM = str(record.CHROM)
    #     BEGIN_POS = str(record.POS)
    #
    #
    #     #If ALT is <INS> or something else --> Make ALT same as POS (Only location of insert needs to be taken into account, not the length)
    #     if  re.search("\w*[\[\]][\w\d]+:\d+[\[\]]\w*" ,str(record.ALT[0])) is not None:
    #         END_CHROM = re.findall("[\w\d]+(?=:)", str(record.ALT[0]))[0]
    #         END_POS = re.findall("(?<=:)[\w\d]+", str(record.ALT[0]))[0]
    #         ID=str(record.ID)
    #
    #         #If ALT_chromosome != POS_chromosome --> Put both breakpoints on new line
    #         if END_CHROM != BEGIN_CHROM:
    #             if record.ALT[0].orientation and record.ALT[0].remoteOrientation:
    #                 print (BEGIN_CHROM + "\t" + str(int(BEGIN_POS)-1) + "\t" + BEGIN_POS + "\t" + ID)
    #                 print (END_CHROM + "\t" + str(int(END_POS)-1) + "\t" + END_POS + "\t" + ID)
    #                 continue
    #             elif not record.ALT[0].orientation and record.ALT[0].remoteOrientation:
    #                 print (BEGIN_CHROM + "\t" + BEGIN_POS + "\t" + str(int(BEGIN_POS)+1) + "\t" + ID)
    #                 print (END_CHROM + "\t" + str(int(END_POS)-1) + "\t" + END_POS + "\t" + ID)
    #                 continue
    #             elif record.ALT[0].orientation and not record.ALT[0].remoteOrientation:
    #                 print (BEGIN_CHROM + "\t" + str(int(BEGIN_POS)-1) + "\t" + BEGIN_POS + "\t" + ID)
    #                 print (END_CHROM + "\t" + END_POS + "\t" + str(int(END_POS)+1) + "\t" + ID)
    #                 continue
    #             elif not record.ALT[0].orientation and not record.ALT[0].remoteOrientation:
    #                 print (BEGIN_CHROM + "\t" + BEGIN_POS + "\t" + str(int(BEGIN_POS)+1) + "\t" + ID)
    #                 print (END_CHROM + "\t" + END_POS + "\t" + str(int(END_POS)+1) + "\t" + ID)
    #                 continue
    #         else:
    #             print (BEGIN_CHROM + "\t" + BEGIN_POS + "\t" + END_POS + "\t" + ID)
    #             #END_COORDINATES = END_CHROM + "\t" + END_POS
    #     elif str(record.ALT[0]) == "<INS>":
    #         #print (str(record.ALT[0]))
    #         END_POS = str(record.INFO["END"])
    #         END_COORDINATES = BEGIN_CHROM + "\t" + str(record.INFO["END"])
    #     else:
    #         print ("ERROR: ALT has an unknown format  --->  ", str(record.ALT[0]))
    #
    #     #if record.CHROM != re.split("[\[\]], record.ALT[0].split(":"))
    #     #print (BEGIN_CHROM + "\t" + BEGIN_POS + "\t" + END_POS)
    #
    #
    #
    #     #print(str(record.CHROM) + "\t" + str(record.POS) + "\t" + str(record.ALT[0]))
    #     #print(pyvcf.model._Breakend)
