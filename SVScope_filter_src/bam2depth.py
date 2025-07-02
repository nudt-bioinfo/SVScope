import subprocess
import os
from utilities import mymkdir

bam_name = "HG002.hs37d5.sorted.bam"
output_file = "output.depth.txt"
data_dir = "/mnt/HHD_16T_1/lyl/data/HG002_GRCh37/new_truvari_results/"

cmd = "samtools depth " + "/mnt/HHD_16T_1/lyl/data/HG002_GRCh37/" + bam_name + " > " + data_dir + output_file
print(cmd)
print("==== starting samtools deal ====")
subprocess.call(cmd, shell = True)

# samtools depth sorted_final_merged.bam > output.depth.txt

mymkdir(data_dir + "depth/")

with open(data_dir + output_file, "r") as f:
    for line in f:
        with open(data_dir + "depth/" + line.split("\t")[0], "a+") as subf:
            subf.write(line)
