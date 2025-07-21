
import codecs
import csv
import os
from enum import Enum
import numpy as np
import torch
import torchvision
import pysam

from bed2image import trans2img

hight = 224
resize = torchvision.transforms.Resize([hight, hight])


def mymkdir(mydir):
    if not os.path.exists(mydir):  
        os.makedirs(mydir)

class IdentifyDataset(torch.utils.data.Dataset):
    def __init__(self, path):

        self.insfile_list = os.listdir(path + "dataset/ins")
        self.delfile_list = os.listdir(path + "dataset/del")
        self.nfile_list = os.listdir(path + "dataset/n")
        self.path = path

        self._len = len(self.insfile_list) + \
            len(self.delfile_list) + len(self.nfile_list)

    def __len__(self):
        return self._len

    def __getitem__(self, index):
        if index < len(self.insfile_list):
            return torch.load(self.path + "dataset/ins/" + self.insfile_list[index])
        elif index < len(self.insfile_list) + len(self.delfile_list):
            index -= len(self.insfile_list)
            return torch.load(self.path + "dataset/del/" + self.delfile_list[index])
        else:
            index -= len(self.insfile_list) + len(self.delfile_list)

            return torch.load(self.path + "dataset/n/" + self.nfile_list[index])

def kernel_cigar(read, ref_min, ref_max, ciagr_resize, zoom):
    # Read-pair types
    SV_SIGNAL_RP_TYPE = Enum("SV_SIGNAL_RP_TYPE", 'LR LLRR RL')

    def get_read_pair_type(read, rl_dist_thr=5):
        # TODO: if one of the reads in the pair is ambiguously mapped, can't infer the correct orientation
        if read.is_reverse == read.mate_is_reverse:
            return SV_SIGNAL_RP_TYPE.LLRR
        if ((read.reference_start + rl_dist_thr) < read.next_reference_start and read.is_read2 and read.is_reverse) or \
                (read.reference_start > (
                        read.next_reference_start + rl_dist_thr) and read.is_read1 and not read.is_reverse):
            # read pairs in R2F1 orientation TODO: add support for R1F2
            return SV_SIGNAL_RP_TYPE.RL
        return SV_SIGNAL_RP_TYPE.LR

    # Provide a default value for read_pair_type
    read_pair_type = None

    cigars_img = torch.zeros([7, int((ref_max - ref_min) / zoom)])

    max_terminal = read.reference_start - ref_min

    for operation, length in read.cigar:  # (operation{10 class}, length)
：
        read_pair_type = get_read_pair_type(read)
        if operation == 0: 
            cigars_img[0, int(max_terminal / zoom):int((max_terminal + length) / zoom)] = 255
            max_terminal += length
        elif operation == 2: 
            cigars_img[1, int(max_terminal / zoom):int((max_terminal + length) / zoom)] = 255
            max_terminal += length
        elif operation == 1: 
            cigars_img[2, int((max_terminal - length / 2) / zoom):int((max_terminal + length / 2) / zoom)] = 255
        elif operation == 4: 
            cigars_img[3, int((max_terminal - length / 2) / zoom):int((max_terminal + length / 2) / zoom)] = 255
        elif operation == 3 or operation == 7 or operation == 8:
            max_terminal += length

    if read_pair_type is not None:

        max_terminal = read.reference_start - ref_min
        read_length = read.reference_end - read.reference_start

        if read_pair_type == SV_SIGNAL_RP_TYPE.LR:
            cigars_img[4, int(max_terminal / zoom):int((max_terminal + read_length) / zoom)] = 100
            # print("successfully!")
        elif read_pair_type == SV_SIGNAL_RP_TYPE.LLRR:
            cigars_img[5, int(max_terminal / zoom):int((max_terminal + read_length) / zoom)] = 150
            # print("successfully!!")
        elif read_pair_type == SV_SIGNAL_RP_TYPE.RL:
            cigars_img[6, int(max_terminal / zoom):int((max_terminal + read_length) / zoom)] = 200
            # print("successfully!!!")
    else:
        pass
        # print("Warning: read_pair_type was not determined for this read.")


    return ciagr_resize(cigars_img.unsqueeze(1))



def cigar_new_img_single_optimal(bam_path, chromosome, begin, end, zoom):
    # print("======= cigar_img_single begin =========")
    r_start = []
    r_end = []
    sam_file = pysam.AlignmentFile(bam_path, "rb")

    for read in sam_file.fetch(chromosome, begin, end):

        if read.reference_start is not None and read.reference_end is not None:
            r_start.append(read.reference_start)
            r_end.append(read.reference_end)

    if r_start :
        ref_min = np.min(r_start)
        ref_max = np.max(r_end)
        cigars_img = torch.empty([7, len(r_start), hight])
        ciagr_resize = torchvision.transforms.Resize([1, hight])

        for i, read in enumerate(sam_file.fetch(chromosome, begin, end)):
            cigars_img[:, i:i + 1,
                       :] = kernel_cigar(read, ref_min, ref_max, ciagr_resize, zoom)

        cigars_img = resize(cigars_img)
    else:
        cigars_img = torch.zeros([7, hight, hight]) 

    sam_file.close()
    # print("======= to input image end =========")
    return cigars_img



