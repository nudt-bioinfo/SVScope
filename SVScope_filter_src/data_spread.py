import utilities as ut
import pandas as pd
import random
import numpy as np
import torch
import torch.nn as nn
from pytorch_lightning.loggers import TensorBoardLogger
import os
from net import IDENet
import pytorch_lightning as pl
from pytorch_lightning import seed_everything
from pytorch_lightning.callbacks import ModelCheckpoint
from multiprocessing import Pool, cpu_count
import pysam
import time
import ray
from ray import tune
from ray.tune import CLIReporter
from ray.tune.search import Repeater
from ray.tune.schedulers import ASHAScheduler, PopulationBasedTraining
from ray.tune.search.hyperopt import HyperOptSearch
from ray.tune.integration.pytorch_lightning import TuneReportCallback, \
    TuneReportCheckpointCallback
import list2img
from hyperopt import hp

data_dir = "/your_data_dir/"
bam_path = "/your_bam_path/" + "HG002.hs37d5.sorted.bam"

ins_vcf_filename = data_dir + "insert_result_data.csv.vcf"
del_vcf_filename = data_dir + "delete_result_data.csv.vcf"

all_enforcement_refresh = 0
position_enforcement_refresh = 0
img_enforcement_refresh = 0
sign_enforcement_refresh = 0  # attention
cigar_enforcement_refresh = 0

# get chr list
sam_file = pysam.AlignmentFile(bam_path, "rb")
chr_list = sam_file.references
chr_length = sam_file.lengths
sam_file.close()

hight = 224

if os.path.exists(data_dir + '/all_n_img' + '.pt'):
    pass
else:
    all_ins_cigar_img = torch.empty(0, 7, hight, hight)
    all_del_cigar_img = torch.empty(0, 7, hight, hight)
    all_negative_cigar_img = torch.empty(0, 7, hight, hight)

    for chromosome, chr_len in zip(chr_list, chr_length):
        print("======= deal " + chromosome + " =======")

        print("position start")
        if os.path.exists(data_dir + 'position/' + chromosome + '/negative' + '.pt') and not position_enforcement_refresh:
            print("loading")
            ins_position = torch.load(
                data_dir + 'position/' + chromosome + '/insert' + '.pt')
            del_position = torch.load(
                data_dir + 'position/' + chromosome + '/delete' + '.pt')
            n_position = torch.load(
                data_dir + 'position/' + chromosome + '/negative' + '.pt')
        else:
            pass
        print("cigar start")
        if os.path.exists(data_dir + 'image/' + chromosome + '/negative_cigar_img' + '.pt') and not cigar_enforcement_refresh:
            print("loading")
            ins_cigar_img = torch.load(
                data_dir + 'image/' + chromosome + '/ins_cigar_img' + '.pt')
            del_cigar_img = torch.load(
                data_dir + 'image/' + chromosome + '/del_cigar_img' + '.pt')
            negative_cigar_img = torch.load(
                data_dir + 'image/' + chromosome + '/negative_cigar_img' + '.pt')
        else:
            pass
        print("cigar end")

        all_ins_cigar_img = torch.cat((all_ins_cigar_img, ins_cigar_img), 0)
        all_del_cigar_img = torch.cat((all_del_cigar_img, del_cigar_img), 0)
        all_negative_cigar_img = torch.cat(
            (all_negative_cigar_img, negative_cigar_img), 0)

    torch.save(all_ins_cigar_img, data_dir + '/all_ins_img' + '.pt')
    torch.save(all_del_cigar_img, data_dir + '/all_del_img' + '.pt')
    torch.save(all_negative_cigar_img, data_dir + '/all_n_img' + '.pt')

print("loading data")

all_ins_img = torch.load(data_dir + '/all_ins_img' + '.pt')
all_del_img = torch.load(data_dir + '/all_del_img' + '.pt')
all_n_img = torch.load(data_dir + '/all_n_img' + '.pt')


print("loaded")

length = len(all_ins_img) + len(all_del_img) + len(all_n_img)

ut.mymkdir(data_dir + "dataset/")
data_dir1 = data_dir + "dataset/"

ut.mymkdir(data_dir1 + "ins/")
ut.mymkdir(data_dir1 + "del/")
ut.mymkdir(data_dir1 + "n/")

for index in range(length):
    print(index)
    if index < len(all_ins_img):
        image = all_ins_img[index].clone()
        torch.save([image, 2], data_dir1 + "ins/" + str(index) + ".pt")
    elif index < len(all_ins_img) + len(all_del_img):
        index -= len(all_ins_img)
        image = all_del_img[index].clone()
        torch.save([image, 1], data_dir1 + "del/" + str(index) + ".pt")
    else:
        index -= len(all_ins_img) + len(all_del_img)
        image = all_n_img[index].clone()
        torch.save([image, 0], data_dir1 + "n/" + str(index) + ".pt")
