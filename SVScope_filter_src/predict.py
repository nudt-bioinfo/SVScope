from itertools import cycle

import matplotlib.pyplot as plt
import numpy as np
import pytorch_lightning as pl
from pytorch_lightning import seed_everything
from scipy import interp
from sklearn.metrics import auc, classification_report, roc_curve

import torch
from net import IDENet
import sys


if len(sys.argv) < 2:
    print("usage error!")
    exit()

# os.environ["CUDA_VISIBLE_DEVICES"] = "3"
seed_everything(2022)

# data_dir = "../datasets/NA12878_PacBio_MtSinai/"
root_dir = "./"

config = {
    "lr": 7.1873e-06,
    "batch_size": 118,  # 14,
    "beta1": 0.9,
    "beta2": 0.999,
    'weight_decay': 0.0011615,
    'model': sys.argv[1]
}


model = IDENet.load_from_checkpoint(
    root_dir + "model/" + config['model'] + ".ckpt", path="/your_model_path/", config=config)

trainer = pl.Trainer(gpus=1)

model.eval()
result = trainer.test(model)

output = torch.load("gold_result.pt")

y = torch.empty(0, 3)
y_hat = torch.empty(0, 3).cuda()


for out in output:
    for ii in out['y']:
        if ii == 0:
            y = torch.cat([y, torch.tensor([1, 0, 0]).unsqueeze(0)], 0)
        elif ii == 1:
            y = torch.cat([y, torch.tensor([0, 1, 0]).unsqueeze(0)], 0)
        else:
            y = torch.cat([y, torch.tensor([0, 0, 1]).unsqueeze(0)], 0)
    y_hat = torch.cat([y_hat, out['y_hat']], 0)


y_test = y.cpu().numpy()
y_score = y_hat.cpu().numpy()
n_classes = y.shape[1]

fpr = dict()
tpr = dict()
roc_auc = dict()
for i in range(n_classes):
    fpr[i], tpr[i], _ = roc_curve(y_test[:, i], y_score[:, i])
    roc_auc[i] = auc(fpr[i], tpr[i])

fpr["micro"], tpr["micro"], _ = roc_curve(y_test.ravel(), y_score.ravel())
roc_auc["micro"] = auc(fpr["micro"], tpr["micro"])

all_fpr = np.unique(np.concatenate([fpr[i] for i in range(n_classes)]))
# Then interpolate all ROC curves at this points
mean_tpr = np.zeros_like(all_fpr)
for i in range(n_classes):
    mean_tpr += interp(all_fpr, fpr[i], tpr[i])
# Finally average it and compute AUC
mean_tpr /= n_classes
fpr["macro"] = all_fpr
tpr["macro"] = mean_tpr
roc_auc["macro"] = auc(fpr["macro"], tpr["macro"])

# Plot all ROC curves
lw = 2
plt.figure()
plt.plot(fpr["micro"], tpr["micro"],
         label='micro-average ROC curve (area = {0:0.2f})'
               ''.format(roc_auc["micro"]),
         color='deeppink', linestyle=':', linewidth=4)

plt.plot(fpr["macro"], tpr["macro"],
         label='macro-average ROC curve (area = {0:0.2f})'
               ''.format(roc_auc["macro"]),
         color='navy', linestyle=':', linewidth=4)

colors = cycle(['aqua', 'darkorange', 'cornflowerblue'])
for i, color in zip(range(n_classes), colors):
    plt.plot(fpr[i], tpr[i], color=color, lw=lw,
             label='ROC curve of class {0} (area = {1:0.2f})'
             ''.format(i, roc_auc[i]))

plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
# plt.title('Some extension of Receiver operating characteristic to multi-class')
plt.legend(loc="lower right")
# plt.savefig("resnet50.pdf", dpi=1000, bbox_inches='tight')
plt.show()


print(classification_report(torch.argmax(y.cpu(), dim=1),
      torch.argmax(y_hat.cpu(), dim=1), digits=4))
