import pandas as pd
import re

data_dir = "/your_data_dir/"


def list_save(filename, data):
    file = open(filename,'w')
    file.writelines(data)
    file.close()
    print(filename + "sucessful")

def set_save(filename, data):
    file = open(filename,'w')
    file.writelines([line+'\n' for line in data])
    file.close()
    print(filename + "sucessful")

insert = ["CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\tSample\tNULL\tNO_CHIP_ID\n"]
delete = ["CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\tSample\tNULL\tNO_CHIP_ID\n"]

filename = data_dir + "unique_to_Tier1.vcf"

chr_list = set()

with open(filename, "r") as f:
    lines = f.readlines()
    for data in lines:
        if "#" in data:
            if "contig=<ID=" in data:
                chr_list.add(re.split("=|,", data)[2])
        else:
            if "SVTYPE=DEL" in data:
                delete.append(data)
            elif "SVTYPE=INS" in data:
                insert.append(data)

list_save(filename + "_ins", insert)
list_save(filename + "_del", delete)
set_save(filename + "_chr", chr_list)

insert_result_data = pd.read_csv(filename + "_ins", sep = "\t")

insert_result_data.insert(2,'SPOS',0)
insert_result_data.insert(3,'EPOS',0)
insert_result_data.insert(4,'SVLEN',0)

for index, row in insert_result_data.iterrows():
    print(index)
    #SPOS, EPOS
    s = row["INFO"]
    pos = s.find("REFWIDENED") # "CIPOS=" POS=POS+6
    svlen = s.find("SVLEN")
    svlen = svlen + 6
    svlen = s[svlen:]
    svlen = svlen.split(";")[0]
    svlen = int(svlen)
    startpos = row["POS"]
    if pos != -1:
        pos = pos + 11 # "REFWIDENED="
        s = s[pos:]
        s = s.split(";")[0]
        # s = s.split(",")
        s = s.split(":")
        s = s[1].split("-")
        # start = int(s[0])
        start = int(s[0])
        # end = int(s[1])
        end = int(s[1])
        insert_result_data.loc[index, ["SPOS"]] = 0
        insert_result_data.loc[index, ["EPOS"]] = 1

    # END
    s = row["INFO"]
    pos = s.find("SVLEN")
    if pos == -1:
        pos = s.find("END") + 4 # "END="
        s = s[pos:]
        s = s.split(";")[0]
        s = int(s) - row["POS"]
        insert_result_data.loc[index, ["SVLEN"]] = s
    else:
        pos = pos + 6 # "SVLEN="
        s = s[pos:]
        s = s.split(";")[0]
        s = int(s)
        insert_result_data.loc[index, ["SVLEN"]] = s


insert_result_data.to_csv(data_dir + "insert_result_data.csv.vcf", sep="\t")



delete_result_data = pd.read_csv(filename + "_del", sep = "\t")

delete_result_data.insert(2,'SPOS',0)
delete_result_data.insert(3,'EPOS',0)
delete_result_data.insert(4,'END',0)
delete_result_data.insert(5,'SEND',0)
delete_result_data.insert(6,'EEND',0)

for index, row in delete_result_data.iterrows():
    print(index)
    #SPOS, EPOS
    s = row["INFO"]
    pos = s.find("REFWIDENED") # "CIPOS=" POS=POS+6
    svlen = s.find("SVLEN")
    svlen = svlen + 6
    svlen = s[svlen:]
    svlen = svlen.split(";")[0]
    svlen = int(svlen)
    startpos = row["POS"]
    if pos != -1:
        pos = pos + 11 # "REFWIDENED="
        s = s[pos:]
        s = s.split(";")[0]
        # s = s.split(",")
        s = s.split(":")
        s = s[1].split("-")
        # start = int(s[0])
        start = int(s[0])
        # end = int(s[1])
        end = int(s[1])
        delete_result_data.loc[index, ["SPOS"]] = 0
        delete_result_data.loc[index, ["EPOS"]] = 1

    # END
    s = row["INFO"]
    pos = s.find("SVLEN")
    if pos == -1:
        pos = s.find("SVLEN") + 6 # "END="
        s = s[pos:]
        s = s.split(";")[0]
        s = row["POS"] - int(s)
        delete_result_data.loc[index, ["END"]] = s
    else:
        pos = pos + 6 # "SVLEN="
        s = s[pos:]
        s = s.split(";")[0]
        s = row["POS"] - int(s)
        delete_result_data.loc[index, ["END"]] = s

    #SPOS, EPOS
    s = row["INFO"]
    pos = s.find("REFWIDENED") # "CIPOS=" POS=POS+6
    svlen = s.find("SVLEN")
    svlen = svlen + 6
    svlen = s[svlen:]
    svlen = svlen.split(";")[0]
    svlen = int(svlen)
    startpos = row["POS"]
    if pos != -1:
        pos = pos + 11 # "REFWIDENED="
        s = s[pos:]
        s = s.split(";")[0]
        # s = s.split(",")
        s = s.split(":")
        s = s[1].split("-")
        # start = int(s[0])
        start = int(s[0])
        # end = int(s[1])
        end = int(s[1])
        delete_result_data.loc[index, ["SEND"]] = 0
        delete_result_data.loc[index, ["EEND"]] = end - startpos + svlen

delete_result_data.to_csv(data_dir + "delete_result_data.csv.vcf", sep="\t")

# index creation
# parallel  samtools index ::: *.bam
