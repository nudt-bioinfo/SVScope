import random
from multiprocessing import cpu_count

import numpy as np
import pytorch_lightning as pl
import torch
import torch.nn as nn
import torch.nn.functional as F
import torchvision
from ray import tune
from sklearn.metrics import classification_report
from torch.utils.data import DataLoader, Subset

import utilities as ut

# model_name = "resnet50"
# model_name = "resnet34"
# model_name = "mobilenet_v2"
# "shufflenet_v2_x1_0"
model_name = "resnet50"
# mnasnet1_0
# resnet50


class attention(nn.Module):
    def __init__(self, dim, out_dim):
        super(attention, self).__init__()
        self.Q_K = nn.Sequential(
            nn.Linear(dim, out_dim),
            nn.Sigmoid(),
        )
        self.V = nn.Sequential(
            nn.Linear(dim, out_dim),
        )

    def forward(self, x):
        qk = self.Q_K(x)
        v = self.V(x)
        out = torch.mul(qk, v)
        return out


class resnet_attention_classfication(nn.Module):
    def __init__(self, full_dim):
        super(resnet_attention_classfication, self).__init__()
        dim1 = full_dim[:-1]
        dim2 = full_dim[1:]
        self.layers = nn.ModuleList(
            nn.Sequential(
                attention(k, m),
                nn.Linear(m, m),
                nn.ReLU(inplace=True),
            ) for k, m in zip(dim1, dim2)
        )

        self.res2 = nn.ModuleList(
            nn.Sequential(
                nn.Linear(full_dim[2 * index], full_dim[2 * (index + 1)]),
                nn.ReLU(inplace=True),
            ) for index in range(int((len(full_dim) - 1) / 2))
        )

    def forward(self, x):
        out = x
        for i in range(len(self.layers)):
            if i % 2 == 0:
                x = out
                out = self.layers[i](out)
            else:
                out = self.layers[i](out) + self.res2[int(i / 2)](x)

        return out

# class MLPClassification(nn.Module):
#     def __init__(self, dims):
#         super(MLPClassification, self).__init__()
#         layers = []
#         dim1 = dims[:-1]
#         dim2 = dims[1:]
#         for in_dim, out_dim in zip(dim1, dim2):
#             layers.append(nn.Linear(in_dim, out_dim))
#             layers.append(nn.ReLU(inplace=True))
#         self.mlp = nn.Sequential(*layers)
#
#     def forward(self, x):
#         return self.mlp(x)

class conv2ds_sequential(nn.Module):
    def __init__(self, full_dim):
        super(conv2ds_sequential, self).__init__()
        dim1 = full_dim[:-1]
        dim2 = full_dim[1:]
        self.layers = nn.ModuleList(
            nn.Sequential(
                nn.Conv2d(in_channels=k, out_channels=m, kernel_size=3,
                          stride=1, padding=1),  # (m, 224, 224)
                nn.BatchNorm2d(m),
                nn.ReLU(inplace=True),
            ) for k, m in zip(dim1, dim2)
        )

    def forward(self, x):
        out = x
        for i, layer in enumerate(self.layers):
            out = layer(out)
        return out


class FocalLoss(nn.Module):
    def __init__(self, alpha=0.25, gamma=2.0, use_sigmoid=True):
        super().__init__()
        self.alpha = alpha
        self.gamma = gamma
        self.use_sigmoid = use_sigmoid
        if use_sigmoid:
            self.sigmoid = nn.Sigmoid()

    def forward(self, pred: torch.Tensor, target: torch.Tensor):
        """
        Focal loss
        :param pred: shape=(B,  HW)
        :param label: shape=(B, HW)
        """
        if self.use_sigmoid:
            pred = self.sigmoid(pred)
        pred = pred.view(-1)
        label = target.view(-1)
        pos = torch.nonzero(label > 0).squeeze(1)
        pos_num = max(pos.numel(), 1.0)
        mask = ~(label == -1)
        pred = pred[mask]
        label = label[mask]
        focal_weight = self.alpha * (label - pred).abs().pow(self.gamma) * (label > 0.0).float(
        ) + (1 - self.alpha) * pred.abs().pow(self.gamma) * (label <= 0.0).float()
        loss = F.binary_cross_entropy(
            pred, label, reduction='none') * focal_weight
        return loss.sum()/pos_num


class IDENet(pl.LightningModule):

    def __init__(self, path, config):
        super(IDENet, self).__init__()

        self.lr = config["lr"]
        self.beta1 = config['beta1']
        self.beta2 = config['beta2']

        self.weight_decay = config['weight_decay']
        self.batch_size = config["batch_size"]
        # model_name = config["model"]

        self.path = path

        conv2d_dim = [7, 7, 3, 3] # [4,4,3,3]修改为[7,7,3,3]
        self.conv2ds = conv2ds_sequential(conv2d_dim)

        # self.resnet_model = eval("torchvision.models." + model_name)(pretrained=True) # [224, 224] -> 1000

        self.resnet_model = eval("torchvision.models." + model_name)(pretrained=False)  # [224, 224] -> 1000

        # state_dict = torch.load("/your_model_path//shufflenetv2_x1-5666bf0f80.pth")
        # self.resnet_model.load_state_dict(state_dict)


        full_dim = [1000, 768, 384, 192, 96, 48, 24, 12, 6]  # test
        self.classfication = resnet_attention_classfication(full_dim)

        self.softmax = nn.Sequential(
            nn.Linear(full_dim[-1], 3),
            nn.Softmax(1)
        )

        # self.criterion = nn.CrossEntropyLoss()
        self.criterion = FocalLoss()  

# class IDENet(pl.LightningModule):
#     def __init__(self, path, config):
#         super(IDENet, self).__init__()
#
#         self.lr = config["lr"]
#         self.beta1 = config['beta1']
#         self.beta2 = config['beta2']
#         self.weight_decay = config['weight_decay']
#         self.batch_size = config["batch_size"]
#         self.path = path
#
#         conv2d_dim = [7, 7, 3, 3]
#         self.conv2ds = conv2ds_sequential(conv2d_dim)
#
#
#         model_name = config["model"]
#         self.resnet_model = eval("torchvision.models." + model_name)(pretrained=True)
#         # state_dict = torch.load(
#         #     "/your_model_path//shufflenetv2_x1-5666bf0f80.pth")
#         # self.resnet_model.load_state_dict(state_dict)
#
#
#         # self.resnet_model.fc = nn.Identity()
#
#
#         full_dim = [1000, 768, 384, 192, 96, 48, 24, 12, 6]
#         self.classfication = MLPClassification(full_dim)
#
#         self.softmax = nn.Sequential(
#             nn.Linear(full_dim[-1], 3),
#             nn.Softmax(1)
#         )
#
#         self.criterion = FocalLoss()
#
#
#


    def training_validation_step(self, batch, batch_idx):
        x, y = batch  # x2(length, 12)
        del batch

        x = self.conv2ds(x)

        # x1 = self.conv2ds(x1[:, 0:1, :, :]) # test

        x = self.resnet_model(x)

        y_t = torch.empty(len(y), 3).cuda()
        for i, y_item in enumerate(y):
            if y_item == 0:
                y_t[i] = torch.tensor([1, 0, 0])
            elif y_item == 1:
                y_t[i] = torch.tensor([0, 1, 0])
            else:
                # y[i] = 1
                y_t[i] = torch.tensor([0, 0, 1])

        y_hat = self.classfication(x)

        y_hat = self.softmax(y_hat)
        loss = self.criterion(y_hat, y_t)
        return loss, y, y_hat

    def training_step(self, batch, batch_idx):
        loss, y, y_hat = self.training_validation_step(batch, batch_idx)

        # logs metrics for each training_step,
        # and the average across the epoch, to the progress bar and logger
        self.log('train_loss', loss, on_step=True,
                 on_epoch=True, prog_bar=True, logger=True)

        # set_trace()
        return {'loss': loss, 'y': y, 'y_hat': torch.argmax(y_hat, dim=1)}

    def training_epoch_end(self, output):
        # set_trace()
        y = []
        y_hat = []

        for out in output:
            y.extend(out['y'])
            y_hat.extend(out['y_hat'])

        y = torch.tensor(y).reshape(-1)
        y_hat = torch.tensor(y_hat).reshape(-1)

        metric = classification_report(y, y_hat, output_dict=True)


        self.log('train_mean', metric['accuracy'], on_step=False,
                 on_epoch=True, prog_bar=True, logger=True)

        # self.log('train_macro_f1', metric['macro avg']['f1-score'], on_step=False, on_epoch=True, prog_bar=True, logger=True)
        self.log('train_macro_pre', metric['macro avg']['precision'],
                 on_step=False, on_epoch=True, prog_bar=True, logger=True)
        self.log('train_macro_re', metric['macro avg']['recall'],
                 on_step=False, on_epoch=True, prog_bar=True, logger=True)

        # self.log('train_0_f1', metric['0']['f1-score'], on_step=False, on_epoch=True, prog_bar=True, logger=True)
        self.log('train_0_pre', metric['0']['precision'],
                 on_step=False, on_epoch=True, prog_bar=True, logger=True)
        self.log('train_0_re', metric['0']['recall'],
                 on_step=False, on_epoch=True, prog_bar=True, logger=True)

        # self.log('train_1_f1', metric['1']['f1-score'], on_step=False, on_epoch=True, prog_bar=True, logger=True)
        self.log('train_1_pre', metric['1']['precision'],
                 on_step=False, on_epoch=True, prog_bar=True, logger=True)
        self.log('train_1_re', metric['1']['recall'],
                 on_step=False, on_epoch=True, prog_bar=True, logger=True)

        # self.log('train_2_f1', metric['2']['f1-score'], on_step=False, on_epoch=True, prog_bar=True, logger=True)
        self.log('train_2_pre', metric['2']['precision'],
                 on_step=False, on_epoch=True, prog_bar=True, logger=True)
        self.log('train_2_re', metric['2']['recall'],
                 on_step=False, on_epoch=True, prog_bar=True, logger=True)

    def validation_step(self, batch, batch_idx):
        loss, y, y_hat = self.training_validation_step(batch, batch_idx)

        # logs metrics for each training_step,
        # and the average across the epoch, to the progress bar and logger
        self.log('validation_loss', loss, on_step=True,
                 on_epoch=True, prog_bar=True, logger=True)
        # set_trace()

        return {'y': y, 'y_hat': torch.argmax(y_hat, dim=1)}

    def validation_epoch_end(self, output):
        y = []
        y_hat = []

        for out in output:
            y.extend(out['y'])
            y_hat.extend(out['y_hat'])

        y = torch.tensor(y).reshape(-1)
        y_hat = torch.tensor(y_hat).reshape(-1)

        metric = classification_report(y, y_hat, output_dict=True)

        self.log('validation_mean', metric['accuracy'],
                 on_step=False, on_epoch=True, prog_bar=True, logger=True)

        # self.log('validation_macro_f1', metric['macro avg']['f1-score'], on_step=False, on_epoch=True, prog_bar=True, logger=True)
        self.log('validation_macro_pre', metric['macro avg']['precision'],
                 on_step=False, on_epoch=True, prog_bar=True, logger=True)
        self.log('validation_macro_re', metric['macro avg']['recall'],
                 on_step=False, on_epoch=True, prog_bar=True, logger=True)

        # self.log('train_0_f1', metric['0']['f1-score'], on_step=False, on_epoch=True, prog_bar=True, logger=True)
        self.log('validation_0_pre', metric['0']['precision'],
                 on_step=False, on_epoch=True, prog_bar=True, logger=True)
        self.log('validation_0_re', metric['0']['recall'],
                 on_step=False, on_epoch=True, prog_bar=True, logger=True)

        # self.log('train_1_f1', metric['1']['f1-score'], on_step=False, on_epoch=True, prog_bar=True, logger=True)
        self.log('validation_1_pre', metric['1']['precision'],
                 on_step=False, on_epoch=True, prog_bar=True, logger=True)
        self.log('validation_1_re', metric['1']['recall'],
                 on_step=False, on_epoch=True, prog_bar=True, logger=True)

        # self.log('train_2_f1', metric['2']['f1-score'], on_step=False, on_epoch=True, prog_bar=True, logger=True)
        self.log('validation_2_pre', metric['2']['precision'],
                 on_step=False, on_epoch=True, prog_bar=True, logger=True)
        self.log('validation_2_re', metric['2']['recall'],
                 on_step=False, on_epoch=True, prog_bar=True, logger=True)

        tune.report(validation_mean=metric['accuracy'])

    def test_step(self, batch, batch_idx):
        loss, y, y_hat = self.training_validation_step(batch, batch_idx)

        # logs metrics for each training_step,
        # and the average across the epoch, to the progress bar and logger
        self.log('validation_loss', loss, on_step=True,
                 on_epoch=True, prog_bar=True, logger=True)
        # set_trace()

        return {'y': y, 'y_hat': y_hat}

    def test_epoch_end(self, output):
        torch.save(output, "gold_result.pt")

    def prepare_data(self):
        train_proportion = 0.8
        input_data = ut.IdentifyDataset(self.path)
        dataset_size = len(input_data)
        indices = list(range(dataset_size))
        split = int(np.floor(train_proportion * dataset_size))
        random.seed(10)
        random.shuffle(indices)
        train_indices, test_indices = indices[:split], indices[split:]
        self.train_dataset = Subset(input_data, train_indices)
        self.test_dataset = Subset(input_data, test_indices)

    def train_dataloader(self):
        # sampler=self.wsampler)
        return DataLoader(dataset=self.train_dataset, batch_size=self.batch_size, pin_memory=True, num_workers=int(cpu_count()), shuffle=True)

    def val_dataloader(self):
        return DataLoader(dataset=self.test_dataset, batch_size=self.batch_size, pin_memory=True, num_workers=int(cpu_count()))

    def test_dataloader(self):
        return DataLoader(dataset=self.test_dataset, batch_size=self.batch_size, pin_memory=True, num_workers=int(cpu_count()))

    # @property
    # def automatic_optimization(self):
    #     return False

    def configure_optimizers(self):
        opt_e = torch.optim.Adam(
            self.parameters(), lr=self.lr, betas=(self.beta1, self.beta2), weight_decay=self.weight_decay)
        # opt_d = torch.optim.Adam(
        #     self.line.parameters(), lr=self.lr, betas=(self.beta1, 0.999))

        return [opt_e]
