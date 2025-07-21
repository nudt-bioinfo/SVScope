import subprocess
import os
from utilities import mymkdir

bam_name = "HG002.hs37d5.sorted.bam"
output_file = "output.depth.txt"
data_dir = "/your_data_dir/"

cmd = "samtools depth " + "/your_bam_path/" + bam_name + " > " + data_dir + output_file
print(cmd)
print("==== starting samtools deal ====")
subprocess.call(cmd, shell = True)

# samtools depth sorted_final_merged.bam > output.depth.txt

mymkdir(data_dir + "depth/")

with open(data_dir + output_file, "r") as f:
    for line in f:
        with open(data_dir + "depth/" + line.split("\t")[0], "a+") as subf:
            subf.write(line)
