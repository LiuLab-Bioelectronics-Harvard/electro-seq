#%%

import tensorflow as tf

import numpy as np
# from utils.dataset import load_bioarxiv_dataset,partitions, Datagen
from model import Model_TE_aug_decoders, custom_build
#new import
import seaborn as sns
import matplotlib.pyplot as plt
import scanpy as sc
import anndata as an
sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.set_figure_params(dpi=80, facecolor='white')
sc.logging.print_header()

import pandas as pd
import pylab
import matplotlib
def color_cm(cmap,NUM_COLORS,):
    color = []
    color_idx = 0
    cm = pylab.get_cmap(cmap)
    for i in range(NUM_COLORS):
        color.append(matplotlib.colors.to_hex(cm(1. * i / NUM_COLORS)))  # color will now be an RGBA tuple
    return color
gene_train_4d_df = pd.read_csv('../SupplementaryCardiacSubstrate/results/csv_data_filtered/gene_data_dpc.csv',header=0,index_col=0)
gene_names = gene_train_4d_df.index
rrr_genes = ['ATP2B1',
 'CACNA1D',
 'HCN4',
 'KCNA4',
 'KCND2',
 'KCND3',
 'KCNH2',
 'KCNJ2',
 'KCNK3',
 'KCNK6',
 'KCNQ1',
 'SLC16A1',
 'SLC25A12',
 'SLC2A1',
 'SLC30A6',
 'SLC44A2',
 'SLC8A1',
 'TRPM4',
 'VDAC2',
 'VDAC3',
 'ATP2A2',
 'LUM']
gene_train_4d_df = gene_train_4d_df.loc[rrr_genes,:]
gene_4d_stage = np.array([vl.split('_')[0] for vl in gene_train_4d_df.columns.values])
gene_4d_stage_unq = np.unique(gene_4d_stage)
gene_4d_cl = color_cm('rainbow',len(gene_4d_stage_unq))
stage_cl = dict(zip(gene_4d_stage_unq,gene_4d_cl))
gene_4d_stage_cl = np.vectorize(stage_cl.get)(gene_4d_stage)
gene_train_4d = gene_train_4d_df.values.T
ephys_train_4d = pd.read_csv('../SupplementaryCardiacSubstrate/results/csv_data_filtered/ephys_data_dpc.csv',header=0,index_col=0).values.T


from sklearn.preprocessing import StandardScaler
scaler_T = StandardScaler()
scaler_E = StandardScaler()
XT = gene_train_4d
XE = ephys_train_4d

median_T = np.median(XT.sum(axis=1))
XT = XT/XT.sum(axis=1)[:, np.newaxis]*median_T
XT = np.log1p(XT)
# XT = XT - np.mean(XT, axis=0)
# XT = XT / np.std(XT, axis=0)
scaler_T.fit(XT)
XT = scaler_T.transform(XT)
scaler_E.fit(XE)
XE = scaler_E.transform(XE)
# XE = XE - np.mean(XE, axis=0)
# XE = XE / np.std(XE, axis=0)

T_dim = XT.shape[1]
E_dim = XE.shape[1]



gene_train_hb_df = pd.read_csv('../SupplementaryExperiment/2_Cardiac Mixed sample/results/gene_data_dpc_device18.csv',header=0,index_col=0)
gene_train_hb_df.columns = gene_names
gene_train_hb_df = gene_train_hb_df[rrr_genes]

gene_hb_stage = np.array(gene_train_hb_df.shape[0]*['day64'])
gene_hb_stage_unq = np.unique(gene_hb_stage)
gene_hb_cl = color_cm('summer',len(gene_hb_stage_unq))
stage_cl = dict(zip(gene_hb_stage_unq,gene_hb_cl))
gene_hb_stage_cl = np.vectorize(stage_cl.get)(gene_hb_stage)

gene_val_hb = gene_train_hb_df.values.T
ephys_val_hb = pd.read_csv('../SupplementaryExperiment/2_Cardiac Mixed sample/results/ephys_data_dpc_device18.csv',header=0,index_col=0).values.T

XT_test = gene_val_hb.T
XE_test = ephys_val_hb.T

XT_test = XT_test/XT_test.sum(axis=1)[:, np.newaxis]*median_T
XT_test = np.log1p(XT_test)
# XT_test = XT_test - np.mean(XT_test, axis=0)
# XT_test = XT_test / np.std(XT_test, axis=0)
XT_test = scaler_T.transform(XT_test)

XE_test = scaler_E.transform(XE_test)
# XE_test = XE_test - np.mean(XE_test, axis=0)
# XE_test = XE_test / np.std(XE_test, axis=0)

#Training step
@tf.function
def train_fn(model, optimizer, xt, xe):
    with tf.GradientTape() as tape:
        zT, zE, XrT, XrE, XrT_cm, XrE_cm = model((xt, xe),
                                train_T=True,
                                train_E=True,
                                augment_decoders=True)

        trainable_weights = [weight for weight in model.trainable_weights]
        loss = sum(model.losses)

    grads = tape.gradient(loss, trainable_weights)
    optimizer.apply_gradients(zip(grads, trainable_weights))
    return zT, zE, XrT, XrE
#Print losses calculated. These are MSE calculations that do not include the
def report_metrics(model, epoch, losstype):
    print('{:10s} Epoch:{:5d}, '
            'mse_T: {:0.3f}, '
            'mse_E: {:0.3f}, '
            'mse_TE: {:0.3f}'.format(losstype,epoch,
                                    model.mse_loss_T.numpy(),
                                    model.mse_loss_E.numpy(),
                                    model.mse_loss_TE.numpy()))
    return model.mse_loss_T.numpy()+model.mse_loss_E.numpy()+model.mse_loss_TE.numpy()
def cal_correlation(XE, XT, XrE, XrT, XrE_cm, XrT_cm):
    ErE = []
    TrT = []
    ErE_cm = []
    TrT_cm = []
    XE = tf.where(tf.math.is_nan(XE), x=0.0, y=XE)
    for i in range(0, XE.shape[0]):
        ErE_cor = np.corrcoef(XE[i, :], XrE[i, :])
        TrT_cor = np.corrcoef(XT[i, :], XrT[i, :])
        ErE_cm_cor = np.corrcoef(XE[i, :], XrE_cm[i, :])
        TrT_cm_cor = np.corrcoef(XT[i, :], XrT_cm[i, :])
        ErE.append(ErE_cor[0, 1])
        TrT.append(TrT_cor[0, 1])
        ErE_cm.append(ErE_cm_cor[0, 1])
        TrT_cm.append(TrT_cm_cor[0, 1])
    print(sum(ErE) / len(ErE), len(ErE))
    print(sum(TrT) / len(TrT), len(TrT))
    print(sum(ErE_cm) / len(ErE_cm), len(ErE_cm))
    print(sum(TrT_cm) / len(TrT_cm), len(TrT_cm))
    cor_ErE = sum(ErE) / len(ErE)
    cor_TrT = sum(TrT) / len(TrT)
    cor_ErE_cm = sum(ErE_cm) / len(ErE_cm)
    cor_TrT_cm = sum(TrT_cm) / len(TrT_cm)
    return cor_ErE,cor_TrT,cor_ErE_cm,cor_TrT_cm


# tf.config.run_functions_eagerly(True)
from utils.dataset import load_bioarxiv_dataset, partitions, Datagen

seed = 42
tf.random.set_seed(seed)
batchsize = 10
n_epochs = 100
n_steps_per_epoch = 1000
maxsteps = n_epochs * n_steps_per_epoch

# Load and build model
cplAE = Model_TE_aug_decoders(T_dim, E_dim, latent_dim=2)
cplAE = custom_build(cplAE, (T_dim, E_dim))

adam_optimizer = tf.keras.optimizers.Adam(learning_rate=1e-4)  # 3)
train_generator = tf.data.Dataset.from_generator(Datagen, output_types=(tf.float32, tf.float32),
                                                 args=(maxsteps, batchsize, XT, XE))
epoch = 0
loss_train = []
loss_val = []
for step, (xT, xE) in enumerate(train_generator):
    train_fn(model=cplAE, optimizer=adam_optimizer, xt=xT, xe=xE)

    # Report training loss
    if (step + 1) % n_steps_per_epoch == 0:
        cplAE((tf.constant(XT, dtype=tf.float32), tf.constant(XE, dtype=tf.float32)),
              train_T=False, train_E=False)
        loss_train.append(report_metrics(model=cplAE, epoch=epoch, losstype='Train'))

        cplAE((tf.constant(XT_test, dtype=tf.float32), tf.constant(XE_test, dtype=tf.float32)),
              train_T=False, train_E=False)
        loss_val.append(report_metrics(model=cplAE, epoch=epoch, losstype='Val.'))

        epoch = epoch + 1






