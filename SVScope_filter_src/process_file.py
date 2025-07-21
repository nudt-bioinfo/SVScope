import argparse
import os
import random

import numpy as np
import pandas as pd
import pytorch_lightning as pl
import torch
import torchvision
from pytorch_lightning import seed_everything
from pytorch_lightning.callbacks import ModelCheckpoint
from pytorch_lightning.loggers import TensorBoardLogger

import utilities as ut

# os.environ["CUDA_VISIBLE_DEVICES"] = "-1"
torch.multiprocessing.set_sharing_strategy('file_system')

seed_everything(2022)

data_dir = "/your_data_dir/"

bam_path = "/your_bam_path/" + "HG002.hs37d5.sorted.bam"

ins_vcf_filename = data_dir + "insert_result_data.csv.vcf"
del_vcf_filename = data_dir + "delete_result_data.csv.vcf"
goldset_ins_vcf_filename = data_dir + "goldset_insert_result_data.csv.vcf"
goldset_del_vcf_filename = data_dir + "goldset_delete_result_data.csv.vcf"

bed_filename = data_dir + "HG002_SVs_Tier1_v0.6.bed"

hight = 224
resize = torchvision.transforms.Resize([hight, hight])

def read_bed_file(bed_filename):
    bed_regions = []
    with open(bed_filename, 'r') as bed_file:
        for line in bed_file:
            fields = line.strip().split()
            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            bed_regions.append((chrom, start, end))
    return bed_regions

bed_regions = read_bed_file(bed_filename)


def position(sum_data):
    chromosome, chr_len = sum_data
    ins_position = []
    del_position = []
    n_position = []

    # insert
    insert_result_data = pd.read_csv(ins_vcf_filename, sep="\t", index_col=0)
    insert_chromosome = insert_result_data[insert_result_data["CHROM"] == chromosome]
    row_pos = []
    for index, row in insert_chromosome.iterrows():
        row_pos.append(row["POS"])


    goldset_insert_result_data = pd.read_csv(goldset_ins_vcf_filename, sep="\t", index_col=0)
    goldset_insert_chromosome = goldset_insert_result_data[goldset_insert_result_data["CHROM"] == chromosome]
    row_pos0 = []
    for index, row in goldset_insert_chromosome.iterrows():
        row_pos0.append(row["POS"])
    set_pos = set()
    for pos in row_pos0:
        set_pos.update(range(pos - 256, pos + 256))

    for pos in row_pos0:
        gap = 1024
        # positive
        begin = pos - 1 - gap
        end = pos - 1 + gap
        if begin < 0:
            begin = 0
        if end >= chr_len:
            end = chr_len - 1

        ins_position.append([begin, end])


    # delete
    delete_result_data = pd.read_csv(del_vcf_filename, sep="\t", index_col=0)
    delete_chromosome = delete_result_data[delete_result_data["CHROM"] == chromosome]
    row_pos = []
    row_end = []
    for index, row in delete_chromosome.iterrows():
        row_pos.append(row["POS"])
        row_end.append(row["END"])

    goldset_delete_result_data = pd.read_csv(goldset_del_vcf_filename, sep="\t", index_col=0)
    goldset_delete_chromosome = goldset_delete_result_data[goldset_delete_result_data["CHROM"] == chromosome]
    row_pos1 = []
    row_end1 = []
    for index, row in goldset_delete_chromosome.iterrows():
        row_pos1.append(row["POS"])
        row_end1.append(row["END"])
    for pos in row_pos1:
        set_pos.update(range(pos - 256, pos + 256))

    for pos, end in zip(row_pos1, row_end1):
        gap = int((end - pos) / 4)
        if gap == 0:
            gap = 1
        # positive
        begin = pos - 1 - gap
        end = end - 1 + gap
        if begin < 0:
            begin = 0
        if end >= chr_len:
            end = chr_len - 1

        del_position.append([begin, end])

        if del_position:
            max_span = max(end - begin for begin, end in del_position)
            print("goldset_del_position max length is:", max_span)
        else:
            print("goldset_del_position is null")

        # 生成 negative
        del_length = end - begin  

        for _ in range(2):
            end = begin

            while end - begin < del_length / 2 + 1:

                bed_region = random.choice(bed_regions)
                if bed_region[0] == chromosome:
                    random_begin = random.randint(bed_region[1], bed_region[2])
                    while random_begin in set_pos:
                        random_begin = random.randint(bed_region[1], bed_region[2])
                    begin = random_begin - 1 - gap
                    end = begin + del_length
                    if begin < 0:
                        begin = 0
                    if end >= chr_len:
                        end = chr_len - 1

            n_position.append([begin, end])

    save_path = data_dir + 'position/' + chromosome

    ut.mymkdir(save_path)
    torch.save(ins_position, save_path + '/insert' + '.pt')
    torch.save(del_position, save_path + '/delete' + '.pt')
    torch.save(n_position, save_path + '/negative' + '.pt')


def create_image(sum_data):
    chromosome, chr_len = sum_data

    print("deal " + chromosome)

    ins_position = torch.load(
        data_dir + 'position/' + chromosome + '/insert' + '.pt')
    del_position = torch.load(
        data_dir + 'position/' + chromosome + '/delete' + '.pt')
    n_position = torch.load(
        data_dir + 'position/' + chromosome + '/negative' + '.pt')

    print("cigar start")

    save_path = data_dir + 'image/' + chromosome

    if os.path.exists(save_path + '/negative_cigar_img' + '.pt'):
        return

    ins_cigar_img = torch.empty(len(ins_position), 7, hight, hight)
    del_cigar_img = torch.empty(len(del_position) * 1, 7, hight, hight)
    negative_cigar_img = torch.empty(len(n_position), 7, hight, hight)

    for i, b_e in enumerate(ins_position):
        # f positive_cigar_img = torch.cat((positive_cigar_img, ut.cigar_img(chromosome_cigar, chromosome_cigar_len, refer_q_table[begin], refer_q_table[end]).unsqueeze(0)), 0)
        zoom = 1
        fail = 1
        while fail:
            try:
                fail = 0
                ins_cigar_img[i] = ut.cigar_new_img_single_optimal(
                    bam_path, chromosome, b_e[0], b_e[1], zoom)
            except Exception as e:
                fail = 1
                zoom += 1
                print(e)
                print("Exception cigar_img_single_optimal " + str(zoom))
                #新增异常捕获机制
                print(f"Exception cigar_img_single_optimal {zoom}: {e}")
                print("Traceback info:")
                print(traceback.format_exc())
                #到此为止

        print("===== finish(ins_cigar_img) " + chromosome + " " + str(i))

    for i, b_e in enumerate(del_position):
        # f positive_cigar_img = torch.cat((positive_cigar_img, ut.cigar_img(chromosome_cigar, chromosome_cigar_len, refer_q_table[begin], refer_q_table[end]).unsqueeze(0)), 0)
        zoom = 1
        fail = 1
        num_images_per_position = 1
        for j in range(num_images_per_position):
            while fail:
                try:
                    fail = 0
                    del_cigar_img[i * num_images_per_position + j] = ut.cigar_new_img_single_optimal(
                        bam_path, chromosome, b_e[0], b_e[1], zoom)
                except Exception as e:
                    fail = 1
                    zoom += 1
                    print(e)
                    print("Exception cigar_img_single_optimal " + str(zoom))

            print("===== finish(del_position) " + chromosome + " " + str(i))

    for i, b_e in enumerate(n_position):
        # f negative_cigar_img = torch.cat((negative_cigar_img, ut.cigar_img(chromosome_cigar, chromosome_cigar_len, refer_q_table[begin], refer_q_table[end]).unsqueeze(0)), 0)
        zoom = 1
        fail = 1
        while fail:
            try:
                fail = 0
                negative_cigar_img[i] = ut.cigar_new_img_single_optimal(
                    bam_path, chromosome, b_e[0], b_e[1], zoom)
            except Exception as e:
                fail = 1
                zoom += 1
                print(e)
                print("Exception cigar_img_single_optimal " + str(zoom))

        print("===== finish(n_position) " + chromosome + " " + str(i))
    ut.mymkdir(save_path) # save_path = data_dir + 'image/' + chromosome
    torch.save(ins_cigar_img, save_path + '/ins_cigar_img' + '.pt')
    torch.save(del_cigar_img, save_path + '/del_cigar_img' + '.pt')
    torch.save(negative_cigar_img, save_path +
               '/negative_cigar_img' + '.pt')
    print("cigar end")


def parse_args():

    description = "you should add those parameter"
    parser = argparse.ArgumentParser(description=description)
    help = "The path of address"
    parser.add_argument('--chr', help=help)
    parser.add_argument('--len', help=help)
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_args()
    # print(args.chr)
    # print(type(args.chr))
    position([args.chr, int(args.len)])
    create_image([args.chr, int(args.len)])
