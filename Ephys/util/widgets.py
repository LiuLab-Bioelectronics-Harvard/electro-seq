import matplotlib.pyplot as plt
from numpy import arange
import numpy as np
from brpylib import NsxFile, brpylib_ver
import pickle
from collections import Counter
import seaborn as sns
from matplotlib.patches import Patch
def read_ns4datafile(datafile):
    brpylib_ver_req = "1.3.1"
    if brpylib_ver.split('.') < brpylib_ver_req.split('.'):
        raise Exception("requires brpylib " + brpylib_ver_req + " or higher, please use latest version")
    elec_ids = 'all'  # 'all' is default for all (1-indexed)                                        
    nsx_file = NsxFile(datafile)
    # Extract data - note: data will be returned based on *SORTED* elec_ids, see cont_data['elec_ids']
    cont_data = nsx_file.getdata(elec_ids)
    unit = nsx_file.extended_headers[1]['Units']
    # Close the nsx file now that all data is out
    nsx_file.close()
    return cont_data, unit
def contingency(a, b, unique_a, unique_b):
    """Populate contingency matrix. Rows and columns are not normalized in any way.
    
    Args:
        a (np.array): labels
        b (np.array): labels
        unique_a (np.array): unique list of labels. Can have more entries than np.unique(a)
        unique_b (np.array): unique list of labels. Can have more entries than np.unique(b)

    Returns:
        C (np.array): contingency matrix.
    """
    assert a.shape == b.shape
    C = np.zeros((np.size(unique_a), np.size(unique_b)))
    for i, la in enumerate(unique_a):
        for j, lb in enumerate(unique_b):
            C[i, j] = np.sum(np.logical_and(a == la, b == lb))
    return C
def grouped_obs_mean(adata, ephys,group_key, cell_integrated,markers,days,clusters,sort=False, layer=None):
    if layer is not None:
        getX = lambda x: x.raw.layers[layer]
    else:
        getX = lambda x: x.raw.X
#     print(adata.obs)
    grouped = adata.obs.groupby(group_key)
    out = pd.DataFrame(
        np.zeros((adata.shape[1], len(grouped)), dtype=np.float64),
        columns=[xx[0]+'_postion_'+str(xx[1]) for xx in list(grouped.groups.keys())],
        index=adata.var_names
    )

    for group, idx in grouped.indices.items():
#         print(group, idx)
        X = getX(adata[idx])
        group_name=group[0]+'_postion_'+str(group[1])
#         out[group_name] = np.log(np.ravel(X.mean(axis=0, dtype=np.float64))+1)
        out[group_name] = np.ravel(X.mean(axis=0, dtype=np.float64))
    out = out.transpose()
    if sort:
        out_idx = out.loc[:, markers].mean().sort_values().index
    else:
        out_idx = markers
    ephys_samples = ephys.index.to_numpy()
    gene_samples = out.index.to_numpy()
    ephys_order=[]
    gene_order=[]
    stage_label_gene=[]
    stage_label_ephy=[]
    t_type_save = []
    e_type_save = []
    if days:
        for stage in days:
            cell_ = cell_integrated[cell_integrated['stage'] == stage].loc[:, ['channel', 'position','cell type','e type (4)']]
            for index, sample in cell_.iterrows():
                ephys_name = stage+'_Channel_'+str(sample['channel'])
                ephys_idx =np.where(ephys_samples == ephys_name)[0]
                if ephys_idx.shape[0]!=0:
                    ephys_order.append(ephys_idx[0])
                    stage_label_ephy.append(stage)
                    e_type_save.append(sample['e type (4)'])
                gene_name = stage+'_postion_'+str(sample['position'])
                
                gene_idx =np.where(gene_samples == gene_name)[0]
                if gene_idx.shape[0]!=0:
                    gene_order.append(gene_idx[0])
                    stage_label_gene.append(stage)
#                     print(gene_idx)
                    t_type_save.append(sample['cell type'])
        network_pal = sns.light_palette('gray', len(np.unique(np.array(stage_label_gene))))
        network_lut = dict(zip(np.unique(np.array(stage_label_gene)), network_pal))
        network_colors_T = tuple(map(tuple, np.array(np.vectorize(network_lut.get)(np.array(stage_label_gene))).T))
        network_colors_E = network_colors_T
        print(len(ephys_order),len(gene_order))
        out_final = out.loc[:, out_idx].iloc[gene_order,:]
        out_save = out.iloc[gene_order,:]
        
        g = sns.clustermap(out_final.T, xticklabels=False,yticklabels=True,col_colors=network_colors_T, row_cluster=False, col_cluster=False,cmap='Reds',standard_scale=0)
        
        handles = [Patch(facecolor=network_lut[name]) for name in network_lut]
        plt.legend(handles, network_lut, title='stage',
        bbox_to_anchor=(1, 1), bbox_transform=plt.gcf().transFigure, loc='upper right')

        plt.show()
        

        ephys_final = ephys.iloc[ephys_order,:]
        x = ephys_final.values #returns a numpy array
        g=sns.clustermap(ephys_final.T, xticklabels=False,yticklabels=False, row_cluster=False,col_cluster=False,cmap='Blues',standard_scale=0)

        plt.show()

    elif clusters:
        print(np.unique(cell_integrated['e type (4)'].values))
        for ctype in np.unique(cell_integrated['cell type'].values):
            cell_t = cell_integrated[cell_integrated['cell type'] == ctype].loc[:, ['channel', 'position','stage','cell type','e type (4)']]
           
            for index, sample in cell_t.iterrows():
                gene_name = sample['stage']+'_postion_'+str(sample['position'])
                gene_idx =np.where(gene_samples == gene_name)[0]
                if gene_idx.shape[0]!=0:
                    gene_order.append(gene_idx[0])
                    stage_label_gene.append(ctype)
                ephys_name = sample['stage']+'_Channel_'+str(sample['channel'])
                ephys_idx =np.where(ephys_samples == ephys_name)[0]
                if ephys_idx.shape[0]!=0:
                    ephys_order.append(ephys_idx[0])
                    stage_label_ephy.append(sample['e type (4)'])
                    
        network_pal_t = sns.light_palette('Gray', len(np.unique(np.array(stage_label_gene))))
        network_pal_e = sns.light_palette('Gray', len(np.unique(np.array(stage_label_ephy))))

        network_lut_t = dict(zip(np.unique(np.array(stage_label_gene)), network_pal_t))
        network_lut_e = dict(zip(np.unique(np.array(stage_label_ephy)), network_pal_e))

        network_colors_T = tuple(map(tuple, np.array(np.vectorize(network_lut_t.get)(np.array(stage_label_gene))).T))
        network_colors_E = tuple(map(tuple, np.array(np.vectorize(network_lut_e.get)(np.array(stage_label_ephy))).T))
        print(out_idx)
        out_final = out.loc[:, out_idx].iloc[gene_order,:]
        out_save = out.iloc[gene_order,:]
        
        ephys_final = ephys.iloc[ephys_order,:]
        g = sns.clustermap(out_final, xticklabels=True,yticklabels=False,row_colors=network_colors_T, row_cluster=False, col_cluster=False,cmap='Reds',standard_scale=1)
        
        handles = [Patch(facecolor=network_lut_t[name]) for name in network_lut_t]
        plt.legend(handles, network_lut_t, title='gene cluster',
        bbox_to_anchor=(1, 1), bbox_transform=plt.gcf().transFigure, loc='upper right')
        plt.show()
        
        g=sns.clustermap(ephys_final, xticklabels=True,yticklabels=False, row_cluster=False,col_cluster=False,cmap='Blues',standard_scale=1)
        plt.show()
        
        
        ephys_order=[]
        gene_order=[]
        stage_label_gene=[]
        stage_label_ephy=[]
        t_type_save = []
        e_type_save = []
        
        for ctype in np.unique(cell_integrated['e type (4)'].values):
            cell_e = cell_integrated[cell_integrated['e type (4)'] == ctype].loc[:, ['channel', 'position','stage','cell type','e type (4)']]
            
            for index, sample in cell_e.iterrows():
                ephys_name = sample['stage']+'_Channel_'+str(sample['channel'])
                ephys_idx =np.where(ephys_samples == ephys_name)[0]
                if ephys_idx.shape[0]!=0:
                    ephys_order.append(ephys_idx[0])
                    stage_label_ephy.append(ctype)
                gene_name = sample['stage']+'_postion_'+str(sample['position'])
                gene_idx =np.where(gene_samples == gene_name)[0]
                if gene_idx.shape[0]!=0:
                    gene_order.append(gene_idx[0])
                    stage_label_gene.append(sample['cell type'])
                    
        network_pal_t = sns.light_palette('Gray', len(np.unique(np.array(stage_label_gene))))
        network_pal_e = sns.light_palette('Gray', len(np.unique(np.array(stage_label_ephy))))

        network_lut_t = dict(zip(np.unique(np.array(stage_label_gene)), network_pal_t))
        network_lut_e = dict(zip(np.unique(np.array(stage_label_ephy)), network_pal_e))

        network_colors_T = tuple(map(tuple, np.array(np.vectorize(network_lut_t.get)(np.array(stage_label_gene))).T))
        network_colors_E = tuple(map(tuple, np.array(np.vectorize(network_lut_e.get)(np.array(stage_label_ephy))).T))
        
        out_final = out.loc[:, out_idx].iloc[gene_order,:]
        out_save = out.iloc[gene_order,:]
        
        ephys_final = ephys.iloc[ephys_order,:]
        g = sns.clustermap(out_final, xticklabels=True,yticklabels=False, row_cluster=False, col_cluster=False,cmap='Reds',standard_scale=1)

        plt.show()
        
        g=sns.clustermap(ephys_final, xticklabels=True,yticklabels=False, row_colors=network_colors_E, row_cluster=False,col_cluster=False,cmap='Blues',standard_scale=1)
        handles = [Patch(facecolor=network_lut_e[name]) for name in network_lut_e]
        plt.legend(handles, network_lut_e, title='ephys cluster',
        bbox_to_anchor=(1, 1), bbox_transform=plt.gcf().transFigure, loc='upper right')
        plt.show()


    print(out_final.shape)
    print(ephys_final.shape)
    return ephys_order,out,ephys_final,out_save,t_type_save,e_type_save
def plot_single_channel(cont_data, plot_chan, unit):
    ch_idx = cont_data['elec_ids'].index(plot_chan)
    t = cont_data['start_time_s'] + arange(cont_data['data'].shape[1]) / cont_data['samp_per_s']
    plt.plot(t, cont_data['data'][ch_idx,:])
    #plt.axis([t[0], t[-1], min(cont_data['data'][ch_idx,:]), max(cont_data['data'][ch_idx,:])])
    plt.locator_params(axis='y', nbins=20)
    plt.xlabel('Time (s)')
    plt.ylabel("Output (" + unit + ")")
    plt.title("Channel " + str(plot_chan))

import pandas as pd
def plot_all_channel(cont_data, unit, sub=True, merge=False, xlim=None):
    if xlim:
        t_start = int(xlim[0] * cont_data['samp_per_s'])
        t_end = int(xlim[1] * cont_data['samp_per_s'])
        t = xlim[0] + arange(t_end - t_start) / cont_data['samp_per_s']
        all_channel = pd.DataFrame(data=cont_data['data'][:, t_start:t_end].transpose())
    else:
        t = cont_data['start_time_s'] + arange(cont_data['data'].shape[1]) / cont_data['samp_per_s']
        all_channel = pd.DataFrame(data=cont_data['data'].transpose())
    all_channel["time(s)"] = t
    if sub:
        all_channel.plot(x="time(s)", ylabel="Output (" + unit + ")", subplots=True, sharey=True, figsize=(10, 10),
                         legend=False)
        plt.legend(loc='right')
        if xlim:
            plt.xlim(xlim)
    if merge:
        all_channel.plot(x="time(s)", ylabel="Output (" + unit + ")", figsize=(10, 5), legend=False)
        plt.legend(loc='right')
        if xlim:
            plt.xlim(xlim)

import copy
from scipy.signal import butter, filtfilt, freqz, iirnotch, welch


def butter_lowpass_filter(data, cutoff, fs=10000, order=5):
    b, a = butter(N=order, Wn=cutoff, btype='low', fs=fs)
    y = filtfilt(b, a, data)
    return y


def iir_notch_filter(data, notch, fs=10000, Q=10.0):
    b, a = iirnotch(notch, Q, fs)
    y = filtfilt(b, a, data, padlen=150)
    return y

def butter_bandpass(lowcut, highcut, fs, order=5):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    return b, a


def butter_bandpass_filter(data, lowcut, highcut, fs, order=5):
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    y = filtfilt(b, a, data)
    y = filtfilt(b, a, data)
    return y



def data_filter(raw_data, cutoff=400, notch=60):
    # Filter requirements.
    order = 6
    fs = raw_data['samp_per_s']  # sample rate, Hz
    # Filter the data, and plot both the original and filtered signals.
    data = raw_data['data']
    for i in range(0, data.shape[0]):
        data[i, :] = butter_lowpass_filter(data[i, :], cutoff, fs, order)
        #data[i, :] = iir_notch_filter(data[i, :], notch, fs, Q=30.0)
    raw_data['data'] = data
    return raw_data
def data_filter_wild(raw_data, fs = 10000, lowcut=None,highcut=400, notch=60):
    # Filter requirements.
    order = 6
    for i in range(0, raw_data.shape[0]):
        if lowcut and highcut:
            raw_data[i, :] = butter_bandpass_filter(raw_data[i, :], lowcut, highcut, fs, order=order)
        else:
            raw_data[i, :] = butter_lowpass_filter(raw_data[i, :], highcut, fs, order=order)
        if notch:
            raw_data[i, :] = iir_notch_filter(raw_data[i, :], notch, fs, Q=10.0)
    return raw_data
## code modified from SpikeInterface
import scipy.signal as ss


# def detect_and_align_peaks_single_channel(cont_data_filtered, channel, n_std, detect_sign, n_pad, min_diff_samples,
#                                           align=True):
#     channel_idx = cont_data_filtered['elec_ids'].index(channel)
#     trace = np.squeeze(cont_data_filtered['data'][channel_idx, :])
#     if detect_sign == -1:
#         thresh = -n_std * np.median(np.abs(trace) / 0.6745)
#         idx_spikes = np.where(trace < thresh)[0]
#     elif detect_sign == 1:
#         thresh = n_std * np.median(np.abs(trace) / 0.6745)
#         idx_spikes = np.where(trace > thresh)[0]
#     else:
#         thresh = n_std * np.median(np.abs(trace) / 0.6745)
#         idx_spikes = np.where((trace > thresh) | (trace < -thresh))[0]
#     intervals = np.diff(idx_spikes)
#     sp_times = []
#
#     for i_t, diff in enumerate(intervals):
#         if diff > min_diff_samples or i_t == len(intervals) - 1:
#             idx_spike = idx_spikes[i_t]
#             if align:
#                 if idx_spike - n_pad > 0 and idx_spike + n_pad < len(trace):
#                     spike = trace[idx_spike - n_pad:idx_spike + n_pad]
#                     t_spike = np.arange(idx_spike - n_pad, idx_spike + n_pad)
#                 elif idx_spike - n_pad < 0:
#                     spike = trace[:idx_spike + n_pad]
#                     spike = np.pad(spike, (np.abs(idx_spike - n_pad), 0), 'constant')
#                     t_spike = np.arange(idx_spike + n_pad)
#                     t_spike = np.pad(t_spike, (np.abs(idx_spike - n_pad), 0), 'constant')
#                 elif idx_spike + n_pad > len(trace):
#                     spike = trace[idx_spike - n_pad:]
#                     spike = np.pad(spike, (0, idx_spike + n_pad - len(trace)), 'constant')
#                     t_spike = np.arange(idx_spike - n_pad, len(trace))
#                     t_spike = np.pad(t_spike, (0, idx_spike + n_pad - len(trace)), 'constant')
#                 if detect_sign == -1:
#                     peak_idx = np.argmin(spike)
#                 elif detect_sign == 1:
#                     peak_idx = np.argmax(spike)
#                 else:
#                     peak_idx = np.argmax(np.abs(spike))
#                 min_time_up = t_spike[peak_idx]
#                 sp_times.append(int(min_time_up))
#             else:
#                 sp_times.append(idx_spike)
#
#     labels = [channel] * len(sp_times)
#
#     return sp_times, labels

from scipy.signal import find_peaks
def detect_and_align_peaks_single_channel(cont_data_filtered, channel, n_std, detect_sign, n_pad, min_diff_samples,
                                          align=True,first_spike=False):
    channel_idx = cont_data_filtered['elec_ids'].index(channel)
    trace = np.squeeze(cont_data_filtered['data'][channel_idx, :])
    if detect_sign == -1:
        thresh = -n_std * np.median(np.abs(trace) / 0.6745)
        idx_spikes = np.where(trace < thresh)[0]
    elif detect_sign == 1:
        thresh = n_std * np.median(np.abs(trace) / 0.6745)
        idx_spikes = np.where(trace > thresh)[0]
    else:
        thresh = n_std * np.median(np.abs(trace) / 0.6745)
        idx_spikes = np.where((trace > thresh) | (trace < -thresh))[0]
    intervals = np.diff(idx_spikes)
    sp_times = []

    for i_t, diff in enumerate(intervals):
        if diff > min_diff_samples or i_t == len(intervals) - 1:
            idx_spike = idx_spikes[i_t]
            if first_spike:
                start = max(0, idx_spike - 6000)
                spike_prior = trace[start:idx_spike]
                # plt.plot(np.arange(int(idx_spike - start)), spike_prior)
                # plt.show()
                if detect_sign==1:
                    height = 1.5*np.median(np.abs(spike_prior) / 0.6745)
                    peaks, _  = find_peaks(spike_prior,height=height,distance=n_pad)
                elif detect_sign==-1:
                    height = 1.5*np.median(np.abs(spike_prior) / 0.6745)
                    peaks, _  = find_peaks(-spike_prior,height=height,distance=n_pad)
                if peaks.size != 0:
                    idx_spike = peaks[0]+start
            if align:
                if idx_spike - n_pad > 0 and idx_spike + n_pad < len(trace):
                    spike = trace[idx_spike - n_pad:idx_spike + n_pad]
                    t_spike = np.arange(idx_spike - n_pad, idx_spike + n_pad)
                elif idx_spike - n_pad < 0:
                    spike = trace[:idx_spike + n_pad]
                    spike = np.pad(spike, (np.abs(idx_spike - n_pad), 0), 'constant')
                    t_spike = np.arange(idx_spike + n_pad)
                    t_spike = np.pad(t_spike, (np.abs(idx_spike - n_pad), 0), 'constant')
                elif idx_spike + n_pad > len(trace):
                    spike = trace[idx_spike - n_pad:]
                    spike = np.pad(spike, (0, idx_spike + n_pad - len(trace)), 'constant')
                    t_spike = np.arange(idx_spike - n_pad, len(trace))
                    t_spike = np.pad(t_spike, (0, idx_spike + n_pad - len(trace)), 'constant')
                if detect_sign == -1:
                    peak_idx = np.argmin(spike)
                elif detect_sign == 1:
                    peak_idx = np.argmax(spike)
                else:
                    peak_idx = np.argmax(np.abs(spike))
                min_time_up = t_spike[peak_idx]
                sp_times.append(int(min_time_up))
            else:
                sp_times.append(idx_spike)

    labels = [channel] * len(sp_times)

    return sp_times, labels


## common mode elimination

# %%

def common_reference(cont_data, ref_type='average', ref_groups=None, channel_ids=None, start_frame=None,
                     end_frame=None, verbose=False, dtype='float32', ref_channel=None):
    channel_idxs = np.array([cont_data['elec_ids'].index(ch) for ch in channel_ids])
    if ref_type == 'median':
        if ref_groups is None:
            if verbose:
                print('Common median reference using all channels')
            traces = cont_data['data'][:, start_frame:end_frame]
            traces = traces - np.median(traces, axis=0, keepdims=True)
            return traces[channel_idxs, :].astype(dtype)
        else:
            new_groups = []
            for g in ref_groups:
                new_chans = []
                for chan in g:
                    if chan in cont_data['elec_ids']:
                        new_chans.append(chan)
                new_groups.append(new_chans)
            if verbose:
                print('Common median in groups: ', new_groups)
            traces = np.vstack(np.array([cont_data['data'][split_group, start_frame:end_frame]
                                         - np.median(cont_data['data'][new_groups, start_frame:end_frame],
                                                     axis=0, keepdims=True) for split_group in new_groups]))
            return traces[channel_idxs, :].astype(dtype)
    elif ref_type == 'average':
        if verbose:
            print('Common average reference using all channels')
        if ref_groups is None:
            traces = cont_data['data'][:, start_frame:end_frame]
            traces = traces - np.mean(traces, axis=0, keepdims=True)
            return traces[channel_idxs, :].astype(dtype)
        else:
            new_groups = []
            for g in ref_groups:
                new_chans = []
                for chan in g:
                    if chan in cont_data['elec_ids']:
                        new_chans.append(chan)
                new_groups.append(new_chans)
            if verbose:
                print('Common average in groups: ', new_groups)
            traces = np.vstack(np.array([cont_data['data'][split_group, start_frame:end_frame]
                                         - np.mean(cont_data['data'][new_groups, start_frame:end_frame],
                                                   axis=0, keepdims=True) for split_group in new_groups]))
            return traces[channel_idxs, :].astype(dtype)
    elif ref_type == 'single':
        if ref_groups is None:
            if verbose:
                print('Reference to channel', ref_channel)
            traces = cont_data['data'][channel_ids, start_frame:end_frame] \
                     - cont_data['data'][ref_channel, start_frame:end_frame]
            return traces.astype(dtype)
        else:
            new_groups = []
            for g in ref_groups:
                new_chans = []
                for chan in g:
                    if chan in cont_data['elec_ids']:
                        new_chans.append(chan)
                new_groups.append(new_chans)
            if verbose:
                print('Reference', new_groups, 'to channels', ref_channel)
            traces = np.vstack(np.array([cont_data['data'][split_group, start_frame:end_frame] \
                                         - cont_data['data'][ref, start_frame:end_frame]
                                         for (split_group, ref) in zip(new_groups, ref_channel)]))
            return traces[channel_idxs, :].astype(dtype)


def common_reference_matrix(data, ref_type='average', ref_groups=None, channel_ids=None, start_frame=None,
                     end_frame=None, verbose=True, dtype='float32', ref_channel=None):
    if ref_type == 'median':
        if ref_groups is None:
            if verbose:
                print('Common median reference using all channels')
            traces = data[:, start_frame:end_frame]
            traces = traces - np.median(traces, axis=0, keepdims=True)
            return traces.astype(dtype)
        else:
            new_groups = []
            for g in ref_groups:
                new_chans = []
                for chan in g:
                    new_chans.append(chan)
                new_groups.append(new_chans)
            if verbose:
                print('Common median in groups: ', new_groups)
            traces = np.vstack(np.array([data[split_group, start_frame:end_frame]
                                         - np.median(data[new_groups, start_frame:end_frame],
                                                     axis=0, keepdims=True) for split_group in new_groups]))
            return traces.astype(dtype)
    elif ref_type == 'average':
        if verbose:
            print('Common average reference using all channels')
        if ref_groups is None:
            traces = data[:, start_frame:end_frame]
            traces = traces - np.mean(traces, axis=0, keepdims=True)
            return traces.astype(dtype)
        else:
            new_groups = []
            for g in ref_groups:
                new_chans = []
                for chan in g:
                    new_chans.append(chan)
                new_groups.append(new_chans)
            if verbose:
                print('Common average in groups: ', new_groups)
            traces = np.vstack(np.array([data[split_group, start_frame:end_frame]
                                         - np.mean(data[new_groups, start_frame:end_frame],
                                                   axis=0, keepdims=True) for split_group in new_groups]))
            return traces.astype(dtype)
    elif ref_type == 'single':
        if ref_groups is None:
            if verbose:
                print('Reference to channel', ref_channel)
            traces = data[channel_ids, start_frame:end_frame] \
                     - data[ref_channel, start_frame:end_frame]
            return traces.astype(dtype)
        else:
            new_groups = []
            for g in ref_groups:
                new_chans = []
                for chan in g:
                    new_chans.append(chan)
                new_groups.append(new_chans)
            if verbose:
                print('Reference', new_groups, 'to channels', ref_channel)
            traces = np.vstack(np.array([data[split_group, start_frame:end_frame] \
                                         - data[ref, start_frame:end_frame]
                                         for (split_group, ref) in zip(new_groups, ref_channel)]))
            return traces.astype(dtype)

## trace whitening

# %%

def get_random_data_for_whitening(cont_data, num_chunks=50, chunk_size=1000, seed=0):
    N = cont_data['data_time_s'] * cont_data['samp_per_s']
    random_ints = np.random.RandomState(seed=seed).randint(0, N - chunk_size, size=num_chunks)
    chunk_list = []
    for ff in random_ints:
        chunk = cont_data['data'][:, ff:ff + chunk_size]
        chunk_list.append(chunk)
    return np.concatenate(chunk_list, axis=1)


def compute_whitening_matrix(cont_data, seed):
    data = get_random_data_for_whitening(cont_data, seed=seed)

    # center the data
    data = data - np.mean(data, axis=1, keepdims=True)

    # Original by Jeremy
    AAt = data @ np.transpose(data)
    AAt = AAt / data.shape[1]
    U, S, Ut = np.linalg.svd(AAt, full_matrices=True)
    W = (U @ np.diag(1 / np.sqrt(S))) @ Ut

    # proposed by Alessio
    # AAt = data @ data.T / data.shape[1]
    # D, V = np.linalg.eig(AAt)
    # W = np.dot(np.diag(1.0 / np.sqrt(D + 1e-10)), V)

    return W


def whitening_trace(cont_data, start_frame, end_frame):
    seed = 0
    whitening_matrix = compute_whitening_matrix(cont_data, seed=seed)
    chunk = cont_data['data'][:, start_frame:end_frame]
    chunk = chunk - np.mean(chunk, axis=1, keepdims=True)
    chunk2 = whitening_matrix @ chunk
    return chunk2

# %% md

## bad channel removal

# %%

def bad_channel_remove(cont_data, bad_threshold=2, bad_channel_ids=None, verbose=True):
    if isinstance(bad_channel_ids, (list, np.ndarray)):
        active_channels = []
        for chan in cont_data['elec_ids']:
            if chan not in bad_channel_ids:
                active_channels.append(chan)
        cont_data['elec_ids'] = active_channels
        cont_data['data'] = cont_data['data'][active_channels - 1, :]
    elif bad_channel_ids is None:
        start_frame = int(cont_data['start_time_s'] * cont_data['samp_per_s'])
        end_frame = start_frame + int(cont_data['data_time_s'] * cont_data['samp_per_s'])
        traces = cont_data['data'][:, start_frame:end_frame]
        stds = np.std(traces, axis=1)
        bad_channel_ids = [ch for ch, std in enumerate(stds) if std > bad_threshold * np.median(stds)]
        if verbose:
            print('Automatically removing channels:', bad_channel_ids)
        active_channels = []
        for chan in cont_data['elec_ids']:
            if chan not in bad_channel_ids:
                active_channels.append(chan)
        cont_data['elec_ids'] = active_channels
        active_channels_idx = [ch_idx - 1 for ch_idx in active_channels]
        cont_data['data'] = cont_data['data'][active_channels_idx, :]
    active_channels = cont_data['elec_ids']
    return cont_data, active_channels
# %% md
import math
import pylab


from statistics import median
from scipy.signal import  detrend
import pywt
## find peaks for different stages
def find_peaks_for_stage(datafile, whiten=False,cutoff=10, filt_bf_avg=1000,
                         rfrac=0.5, n_std_min=5, whole_plot=False, average_plot=True,
                         whole_xlim=[0, 60], average_xlim=[0, 0.2], n_pad_left=1000,
                         n_pad_right=9000,n_pad_left_sub=100, n_pad_right_sub=900,cmap='rainbow',
                         n_pad=300,first_spike=False,single_spike=False,save_average_file=False,save_channels=range(1, 16),plot_start=False):
    cont_data, unit = read_ns4datafile(datafile)
    cont_data['data'] = cont_data['data'][:16,:]
    cont_data['elec_ids'] = cont_data['elec_ids'][:16]
    cont_data['ExtendedHeaderIndices'] = cont_data['ExtendedHeaderIndices'][:16]
    cont_data_peak_find = copy.deepcopy(cont_data)
    # dwave = 'dmey'
    # mode = 'symmetric'
    # signal = cont_data_peak_find['data'][0, :]
    # cD = pywt.wavedec(signal, dwave, mode=mode, level=None)
    # test5 = pywt.waverec(cD[0:8], dwave, mode=mode)

    # cont_data_peak_find['data'] = np.squeeze(
    #     common_reference(cont_data_peak_find, channel_ids=cont_data_peak_find['elec_ids'], ref_type='median'))
    cont_data_peak_find = data_filter(cont_data_peak_find, cutoff=cutoff)
    cont_data_peak_find['data'] = np.diff(cont_data_peak_find['data'], axis=1)
    #cont_data_peak_find['data'] = detrend(cont_data_peak_find['data'],type = 'cons', axis = 0)

    # trace whitening
    if whiten:
        start_frame = int(cont_data_peak_find['start_time_s'] * cont_data_peak_find['samp_per_s'])
        end_frame = start_frame + int(cont_data['data_time_s'] * cont_data_peak_find['samp_per_s'])
        cont_data_peak_find['data'] = whitening_trace(cont_data_peak_find, start_frame, end_frame)
    min_dis = rfrac * cont_data['samp_per_s']
    start_peaks = {'time':[],'channel':[]}
    t_spike = arange(n_pad_left + n_pad_right) / cont_data['samp_per_s']
    averaged_spikes = np.zeros([len(cont_data['elec_ids']), n_pad_left + n_pad_right])
    all_spikes = []
    if save_channels:
        NUM_COLORS = len(save_channels)
    else:
        NUM_COLORS = len(cont_data['elec_ids'])
    color = []
    color_idx = 0
    cm = pylab.get_cmap(cmap)
    for i in range(NUM_COLORS):
        color.append(cm(1. * i / NUM_COLORS))  # color will now be an RGBA tuple
    peaks_group = []
    for ch in cont_data['elec_ids']:
        channel_idx = cont_data['elec_ids'].index(ch)
        peak_pos, _ = detect_and_align_peaks_single_channel(cont_data_peak_find, ch, n_std_min,
                                                                1, n_pad, min_dis,first_spike=first_spike)
        peak_neg, _ = detect_and_align_peaks_single_channel(cont_data_peak_find, ch, n_std_min,
                                                                -1, n_pad, min_dis,first_spike=first_spike)
        if (not peak_pos) and (not peak_neg):
            cont_data, _ = bad_channel_remove(cont_data, bad_channel_ids=ch)
            cont_data_peak_find, _ = bad_channel_remove(cont_data_peak_find, bad_channel_ids=ch)
            #print('removed bad channel ' + str(ch))
            bad_spikes = np.zeros([1, n_pad_left + n_pad_right])
            all_spikes.append(bad_spikes)
            peaks_group.append([])
            continue
        if filt_bf_avg:
            cont_data = data_filter(cont_data, cutoff=filt_bf_avg)
        peak_pos_averaged_pos, peak_pos_pos = spike_average(cont_data['data'][channel_idx, :], peak_pos, detect_sign=1,
                                                          n_pad_left=n_pad_left, n_pad_right=n_pad_right,
                                                          n_pad_left_sub=n_pad_left_sub,
                                                          n_pad_right_sub=n_pad_right_sub)
        peak_pos_averaged_neg, peak_pos_neg = spike_average(cont_data['data'][channel_idx, :], peak_pos, detect_sign=-1,
                                                          n_pad_left=n_pad_left, n_pad_right=n_pad_right,
                                                          n_pad_left_sub=n_pad_left_sub,
                                                          n_pad_right_sub=n_pad_right_sub)
        peak_neg_averaged_pos, peak_neg_pos = spike_average(cont_data['data'][channel_idx, :], peak_neg, detect_sign=1,
                                                          n_pad_left=n_pad_left, n_pad_right=n_pad_right,
                                                          n_pad_left_sub=n_pad_left_sub,
                                                          n_pad_right_sub=n_pad_right_sub)
        peak_neg_averaged_neg, peak_neg_neg = spike_average(cont_data['data'][channel_idx, :], peak_neg, detect_sign=-1,
                                                          n_pad_left=n_pad_left, n_pad_right=n_pad_right,
                                                          n_pad_left_sub=n_pad_left_sub,
                                                          n_pad_right_sub=n_pad_right_sub)
        cor_pos_pos = spikes_coef(peak_pos_pos, peak_pos_averaged_pos)
        cor_pos_neg = spikes_coef(peak_pos_neg, peak_pos_averaged_neg)
        if len(peak_pos)<=1:
            cor_pos_pos = cor_pos_pos/2
            cor_pos_neg = cor_pos_neg/2
        cor_neg_pos = spikes_coef(peak_neg_pos, peak_neg_averaged_pos)
        cor_neg_neg = spikes_coef(peak_neg_neg, peak_neg_averaged_neg)
        if len(peak_neg)<=1:
            cor_neg_pos = cor_neg_pos/2
            cor_neg_neg = cor_neg_neg/2
        cor_list = [cor_pos_pos,cor_pos_neg,cor_neg_pos,cor_neg_neg]

        cor_list = [0 if math.isnan(x) else x for x in cor_list]
        #print(cor_list)
        cor_max = cor_list.index(max(cor_list))
        if cor_max == 0:
            averaged_spike = peak_pos_averaged_pos
            spikes = peak_pos_pos
            peak = peak_pos
        elif cor_max == 1:
            averaged_spike = peak_pos_averaged_neg
            spikes = peak_pos_neg
            peak = peak_pos
        elif cor_max == 2:
            averaged_spike = peak_neg_averaged_pos
            spikes = peak_neg_pos
            peak = peak_neg
        else:
            averaged_spike = peak_neg_averaged_neg
            spikes = peak_neg_neg
            peak = peak_neg
        if peak:
            peaks_group.append(peak)
        else:
            peaks_group.append([])
        if single_spike:
            if spikes.shape[0]>1:
                averaged_spike = spikes[1, :]
        averaged_spikes[ch - 1, :] = averaged_spike

        all_spikes.append(spikes)
        if whole_plot:
            plt.figure()
            plot_single_channel(cont_data_peak_find, ch, unit)
            plt.scatter(np.array(peak) / cont_data_peak_find['samp_per_s'],
                        cont_data_peak_find['data'][channel_idx, peak],
                        marker="v", c='red', s=40)
            plt.xlim(whole_xlim)
            plt.show()
            plt.figure()
            plot_single_channel(cont_data, ch, unit)
            plt.scatter(np.array(peak) / cont_data['samp_per_s'], cont_data['data'][channel_idx, peak],
                        marker="v", c='red', s=40)
            plt.xlim(whole_xlim)
            plt.show()
        if average_plot:
            plt.figure(figsize=(0.6,0.4))
            plt.plot(t_spike, spikes.transpose(), color='0.5', alpha=0.2)
            if save_average_file and channel_idx in save_channels:
                plt.plot(t_spike, averaged_spike, linewidth=0.3,color=color[color_idx])
                #plt.title(['channel' + str(ch) + ' averaged'])
                plt.xlim(average_xlim)
                plt.axis('off')
                plt.savefig(save_average_file+'channel'+str(ch)+'.pdf', transparent=True)
                color_idx = color_idx+1
            else:
                plt.plot(t_spike, averaged_spike, linewidth=0.3)
                # plt.title(['channel' + str(ch) + ' averaged'])
                plt.xlim(average_xlim)
                plt.axis('off')
                plt.show()
        if peak[0]==0:
            continue
        #if cont_data['data'][channel_idx, peak[0]]<1:
        loc_start = peak[0]-200
        if loc_start<0:
            loc_start=0
        loc_end = peak[0]+200
        peak_start = np.argmax(cont_data['data'][channel_idx, loc_start:loc_end])+loc_start
        start_peaks['time'].append(peak_start)
        start_peaks['channel'].append(channel_idx)

    start_peaks['channel'] = [x for _, x in sorted(zip(start_peaks['time'], start_peaks['channel']))]
    start_peaks['time'] = sorted(start_peaks['time'])
    if plot_start:
        bad_start1 = np.where((np.diff(start_peaks['time'])) == 0)[0]+1
        bad_start2 = np.where((start_peaks['time']-start_peaks['time'][0])> int(median(start_peaks['time'])))[0]
        good_start = np.setdiff1d(np.arange(len(start_peaks['time'])), np.unique(np.append(bad_start1, bad_start2)))
        start_peaks['time'] = [start_peaks['time'][x] for x in good_start]
        start_peaks['channel'] = [start_peaks['channel'][x] for x in good_start]
        center_time = int(median(start_peaks['time']))
        t_start = center_time - 1000
        if t_start<0:
            t_start=0
        t_end = center_time+1000
        t  = arange(t_end-t_start) / cont_data['samp_per_s']
        all_channel = pd.DataFrame(data=cont_data['data'][start_peaks['channel'],t_start:t_end].transpose())
        all_channel["time(s)"] = t
        all_channel.plot(x="time(s)", ylabel="Output (" + unit + ")", figsize=(10, 5), legend=False,cmap=cmap)
    return averaged_spikes, all_spikes,start_peaks,peaks_group

# %%
def find_peaks_stage2(datafile, whiten=False, detect_sign=1, rfrac=0.5, n_std_min=5, whole_plot=False, average_plot=True,
               whole_xlim=[0, 60], average_xlim=[0, 0.2], n_pad_left=1000, n_pad_right=9000, n_pad_left_sub=100,
               n_pad_right_sub=900):
    cont_data, unit = read_ns4datafile(datafile)
    cont_data_peak_find = copy.deepcopy(cont_data)
    cont_data_peak_find['data'] = np.squeeze(
        common_reference(cont_data_peak_find, channel_ids=cont_data_peak_find['elec_ids'], ref_type='median'))
    cont_data_peak_find = data_filter(cont_data_peak_find, cutoff=100)
    cont_data_peak_find['data'] = np.diff(cont_data_peak_find['data'], axis=1)
    # cont_data_peak_find['data'] = signal.detrend(cont_data_peak_find['data'],type = 'linear', axis = 0)

    # trace whitening
    if whiten:
        start_frame = int(cont_data_peak_find['start_time_s'] * cont_data_peak_find['samp_per_s'])
        end_frame = start_frame + int(cont_data['data_time_s'] * cont_data_peak_find['samp_per_s'])
        cont_data_peak_find['data'] = whitening_trace(cont_data_peak_find, start_frame, end_frame)
    min_dis = rfrac * cont_data['samp_per_s']
    peaks = []
    n_pad = 3000
    detect_sign = detect_sign
    t_spike = arange(n_pad_left + n_pad_right) / cont_data['samp_per_s']
    averaged_spikes = np.zeros([len(cont_data['elec_ids']), n_pad_left + n_pad_right])

    for ch in cont_data['elec_ids']:
        channel_idx = cont_data['elec_ids'].index(ch)
        peak, labels = detect_and_align_peaks_single_channel(cont_data_peak_find, ch, n_std_min,
                                                                detect_sign, n_pad, min_dis)
        if not peak:
            cont_data, _ = bad_channel_remove(cont_data, bad_channel_ids=ch)
            cont_data_peak_find, _ = bad_channel_remove(cont_data_peak_find, bad_channel_ids=ch)
            print('removed bad channel ' + str(ch))
            continue
        peak = [peak_ + 1 for peak_ in peak]

        averaged_spike_pos, spikes_pos = spike_average(cont_data['data'][channel_idx, :], peak, detect_sign=1,
                                                          n_pad_left=n_pad_left, n_pad_right=n_pad_right,
                                                          n_pad_left_sub=n_pad_left_sub,
                                                          n_pad_right_sub=n_pad_right_sub)
        averaged_spike_neg, spikes_neg = spike_average(cont_data['data'][channel_idx, :], peak, detect_sign=-1,
                                                          n_pad_left=n_pad_left, n_pad_right=n_pad_right,
                                                          n_pad_left_sub=n_pad_left_sub,
                                                          n_pad_right_sub=n_pad_right_sub)
        cor_pos = spikes_coef(spikes_pos, averaged_spike_pos)
        cor_neg = spikes_coef(spikes_neg, averaged_spike_neg)
        if cor_pos > cor_neg:
            averaged_spike = averaged_spike_pos
            spikes = spikes_pos
        else:
            averaged_spike = averaged_spike_neg
            spikes = spikes_neg
        averaged_spikes[ch - 1, :] = averaged_spike
        if whole_plot:
            plt.figure()
            plot_single_channel(cont_data_peak_find, ch, unit)
            plt.scatter(np.array(peak) / cont_data_peak_find['samp_per_s'],
                        cont_data_peak_find['data'][channel_idx, peak],
                        marker="v", c='red', s=40)
            plt.xlim(whole_xlim)
            plt.show()
            plt.figure()
            plot_single_channel(cont_data, ch, unit)
            plt.scatter(np.array(peak) / cont_data['samp_per_s'], cont_data['data'][channel_idx, peak],
                        marker="v", c='red', s=40)
            plt.xlim(whole_xlim)
            plt.show()
        if average_plot:
            plt.figure()
            averaged_spike_plot(spikes, averaged_spike, t_spike)
            plt.title(['channel' + str(ch) + ' averaged'])
            plt.xlim(average_xlim)
            plt.show()
        peaks.append(peak)
    return averaged_spikes

def spike_average(trace, sp_times, detect_sign = -1, n_pad_left=1000, n_pad_right=9000,n_pad_left_sub=100, n_pad_right_sub=900):
    elec_num = len(sp_times)
    spikes = np.zeros([elec_num, n_pad_left + n_pad_right])
    for i_t, idx_spike in enumerate(sp_times):
        #padded_ = 0
        if idx_spike - n_pad_left > 0 and idx_spike + n_pad_right < len(trace):
            spike = trace[idx_spike - n_pad_left:idx_spike + n_pad_right]
            t_spike = np.arange(idx_spike - n_pad_left, idx_spike + n_pad_right)
        elif idx_spike - n_pad_left < 0:
            spike = trace[:idx_spike + n_pad_right]
            spike = np.pad(spike, (np.abs(idx_spike - n_pad_left), 0), 'constant')
            t_spike = np.arange(idx_spike + n_pad_right)
            t_spike = np.pad(t_spike, (np.abs(idx_spike - n_pad_left), 0), 'constant')
            #padded_ = np.abs(idx_spike - n_pad_left)
        elif idx_spike + n_pad_right > len(trace):
            spike = trace[idx_spike - n_pad_left:]
            spike = np.pad(spike, (0, idx_spike + n_pad_right - len(trace)), 'constant')
            t_spike = np.arange(idx_spike - n_pad_left, len(trace))
            t_spike = np.pad(t_spike, (0, idx_spike + n_pad_right - len(trace)), 'constant')
        if detect_sign == -1:
            peak_idx = np.argmin(spike[n_pad_left - n_pad_left_sub:n_pad_left + n_pad_right_sub])
        elif detect_sign == 1:
            peak_idx = np.argmax(spike[n_pad_left - n_pad_left_sub:n_pad_left + n_pad_right_sub])
        else:
            peak_idx = np.argmax(np.abs(spike[n_pad_left - n_pad_left_sub:n_pad_left + n_pad_right_sub]))
        peak_idx = idx_spike - (n_pad_left_sub - peak_idx)
        if peak_idx - n_pad_left > 0 and peak_idx + n_pad_right < len(trace):
            spike = trace[peak_idx - n_pad_left:peak_idx + n_pad_right]
            t_spike = np.arange(peak_idx - n_pad_left, peak_idx + n_pad_right)
        elif peak_idx - n_pad_left < 0:
            spike = trace[:peak_idx + n_pad_right]
            spike = np.pad(spike, (np.abs(peak_idx - n_pad_left), 0), 'constant')
            t_spike = np.arange(peak_idx + n_pad_right)
            t_spike = np.pad(t_spike, (np.abs(peak_idx - n_pad_left), 0), 'constant')
        elif peak_idx + n_pad_right > len(trace):
            spike = trace[peak_idx - n_pad_left:]
            spike = np.pad(spike, (0, peak_idx + n_pad_right - len(trace)), 'constant')
            t_spike = np.arange(peak_idx - n_pad_left, len(trace))
            t_spike = np.pad(t_spike, (0, peak_idx + n_pad_right - len(trace)), 'constant')
        spikes[i_t, :] = spike
    averaged_spike = np.mean(spikes, axis=0)
    return averaged_spike, spikes
import os
from features_ECG import *
def extract_all_spikes(data_folder, cutoff,rfrac, filt_bf_avg, n_std_min, whole_plot,
                         average_plot,whole_xlim,average_xlim, n_pad_left,
                         n_pad_right,n_pad_left_sub,n_pad_right_sub, n_pad,
                         first_spike,single_spike,plot_start,cmap,save_average_file=False,spe='002'):
    averaged_spikes = []
    all_spikes = []
    starting_spikes = []
    flags = []
    dir_single_group = []
    stage_dir_name = os.walk(data_folder)
    peak_RR=[]
    for path,_,sub_file_list in stage_dir_name:
        for file in sub_file_list:
            if file.endswith(spe + ".ns4"):
                dir_single_group.append(os.path.join(path,file))
                flags.append(file[-9])
        #print(dir_single_group)
        #print(flags)
        flags_sort_idx = sorted(range(len(flags)), key=lambda k: flags[k])
        #print(flags_sort_idx)
        flags = []
        for idx in flags_sort_idx:
            #print(dir_single_group[idx])
            averaged_spike,spikes,starting_spike,peak = find_peaks_for_stage(dir_single_group[idx], cutoff=cutoff,rfrac=rfrac, filt_bf_avg=filt_bf_avg,
                                         n_std_min=n_std_min, whole_plot=whole_plot,
                                         average_plot=average_plot,whole_xlim=whole_xlim,
                                         average_xlim=average_xlim, n_pad_left=n_pad_left,
                                         n_pad_right=n_pad_right,n_pad_left_sub=n_pad_left_sub,
                                         n_pad_right_sub=n_pad_right_sub,
                                         n_pad=n_pad,cmap=cmap,save_average_file=save_average_file,
                                         first_spike=first_spike,single_spike=single_spike,plot_start=plot_start)
            peak_RR.append([compute_RR_intervals(single_trace) for single_trace in peak])
            all_spikes.append(spikes)
            averaged_spikes.append(averaged_spike)
            starting_spikes.append(starting_spike)
    peak_RR = [item for sublist in peak_RR for item in sublist]
    peak_RR_features = np.array([], dtype=float)

    print("Computing RR intervals ...")

    f_RR = np.empty((0, 4))
    for p in range(len(peak_RR)):
        if peak_RR[p]==[]:
            f_RR=np.vstack((f_RR, np.zeros((1, 4))))
        else:
            avg_pre_R = np.average(peak_RR[p].pre_R)
            avg_post_R = np.average(peak_RR[p].post_R)
            avg_local_R = np.average(peak_RR[p].local_R)
            avg_global_R = np.average(peak_RR[p].global_R)
            f_RR = np.vstack((f_RR, np.array([avg_pre_R,avg_post_R,avg_local_R,avg_global_R])))
    peak_RR_features = np.column_stack((peak_RR_features, f_RR)) if peak_RR_features.size else f_RR
    for i in range(peak_RR_features.shape[1]):
        idx = np.logical_or(peak_RR_features[:, i] > peak_RR_features[:, i].mean(), peak_RR_features[:, i] == 0)
        peak_RR_features[idx, i] = peak_RR_features[:, i].mean()
    averaged_spikes = np.concatenate( averaged_spikes, axis=0 )

    return averaged_spikes,all_spikes,starting_spikes,peak_RR_features

def extract_avraged_spikes(data_folder, cutoff,rfrac, filt_bf_avg, n_std_min, whole_plot,
                         average_plot,whole_xlim,average_xlim, n_pad_left,
                         n_pad_right,n_pad_left_sub,n_pad_right_sub, n_pad,
                         first_spike,single_spike):
    averaged_spikes = []
    flags = []
    dir_single_group = []
    stage_dir_name = os.walk(data_folder)
    for path,_,sub_file_list in stage_dir_name:
        for file in sub_file_list:
            if file.endswith("002.ns4"):
                dir_single_group.append(os.path.join(path,file))
                flags.append(file[-9])
        #print(dir_single_group)
        #print(flags)
        flags_sort_idx = sorted(range(len(flags)), key=lambda k: flags[k])
        #print(flags_sort_idx)
        flags = []
        for idx in flags_sort_idx:
            #print(dir_single_group[idx])
            averaged_spikes.append(find_peaks_for_stage(dir_single_group[idx], cutoff=cutoff,rfrac=rfrac, filt_bf_avg=filt_bf_avg,
                                         n_std_min=n_std_min, whole_plot=whole_plot,
                                         average_plot=average_plot,whole_xlim=whole_xlim,
                                         average_xlim=average_xlim, n_pad_left=n_pad_left,
                                         n_pad_right=n_pad_right,n_pad_left_sub=n_pad_left_sub,
                                         n_pad_right_sub=n_pad_right_sub, n_pad=n_pad,
                                         first_spike=first_spike,single_spike=single_spike))
    #averaged_spikes = np.concatenate( averaged_spikes, axis=0 )
    return averaged_spikes



def align_spikes(spikes,new_samp=3000,detect_sign=1):
    center = new_samp/10
    ch_num = spikes.shape[0]
    time =spikes.shape[1]
    for i in range(ch_num):
        if detect_sign == -1:
            peak_idx = np.argmin(spikes[i,:])
        elif detect_sign == 1:
            peak_idx = np.argmax(spikes[i,:])
        else:
            peak_idx = np.argmax(np.abs(spikes[i,:]))
        diff = int(peak_idx - center)
        if diff < 0:
            spike = spikes[i,:time+diff]
            spike = np.pad(spike, (np.abs(diff), 0), 'constant')
        elif diff > 0:
            spike = spikes[i,diff:]
            spike = np.pad(spike, (0, diff), 'constant')
        else:
            spike = spikes[i,:]
        spikes[i, :] = spike
    return spikes

def averaged_spike_plot(spikes, averaged_spike, t):
    plt.plot(t, spikes.transpose(), color='0.5', alpha=0.2)
    ax = plt.plot(t, averaged_spike, linewidth=2)
    #plt.title(['channel'+ str(channel) +' averaged'])
    return ax
def save_obj(obj, name ):
    with open('obj/'+ name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_obj(name ):
    with open('obj/' + name + '.pkl', 'rb') as f:
        return pickle.load(f)

def vcorrcoef(X,y):
    Xm = np.reshape(np.mean(X,axis=1),(X.shape[0],1))
    ym = np.mean(y)
    r_num = np.sum((X-Xm)*(y-ym),axis=1)
    r_den = np.sqrt(np.sum((X-Xm)**2,axis=1)*np.sum((y-ym)**2))
    r = r_num/r_den
    return r
def spikes_coef(spikes, averaged_spike):
    r = vcorrcoef(spikes, averaged_spike)
    return np.mean(r)

from scipy.signal import resample
def down_sampling(averaged_spikes, sampling_size_new = 3000,ch_num=None):
    start_time = 0
    end_time = averaged_spikes.shape[1]
    t_spike = arange(end_time) / 10000
    averaged_spikes_dsamp = resample(averaged_spikes, sampling_size_new, axis=1)
    print(end_time/ 10000)
    sampling_rate_new = sampling_size_new/(end_time/ 10000)

    t_spike_new = arange(sampling_size_new) / sampling_rate_new
    if ch_num:
        for i in range(ch_num):
            plt.figure()
            plt.plot(t_spike, averaged_spikes[i,start_time:end_time])
            plt.title(['channel'+ str(i) +' averaged before filtering'])
            plt.show()
            #averaged_spikes_filtered[i, :] = resample(averaged_spikes_filtered[i, :], 6000)
            plt.figure()
            plt.plot(t_spike_new, averaged_spikes_dsamp[i,:])
            plt.title(['channel'+ str(i) +' averaged after filtering'])
            plt.show()
    return averaged_spikes_dsamp, sampling_rate_new
#%%

from scipy.signal import welch
import pandas as pd
def cm_fre(averaged_spikes):
    t_spike = arange(10000) / 10000
    average_xlim=[0,1]
    all_channel = pd.DataFrame(data=averaged_spikes.transpose())
    all_channel["time(s)"]=t_spike
    all_channel.plot(x="time(s)",ylabel = "Output (uV)",figsize=(10,5), legend=False)
    plt.title('channels averaged')
    plt.xlim(average_xlim)
    plt.show()
    start_time = 800
    end_time = 1300
    t_spike = arange(start_time,end_time) / 10000
    averaged_spikes_filtered = copy.deepcopy(averaged_spikes)
    for i in range(10):
        plt.figure()
        f_raw, p_raw = welch(averaged_spikes[i,start_time:end_time], fs=10000 )
        plt.semilogy(f_raw, p_raw, label='before filtering')
        plt.xlabel('Frequency (Hz)')
        plt.ylabel('spectral energy density')
        print(f_raw[np.argmax(p_raw)])
        #print(f_raw[np.where(p_raw==p_raw.max())],p_raw.max())
        plt.figure()
        plt.plot(t_spike, averaged_spikes[i,start_time:end_time])
        plt.title(['channel'+ str(i) +' averaged before filtering'])
        plt.show()
        averaged_spikes_filtered[i, :] = iir_notch_filter(averaged_spikes[i, :], 39.0625, fs=10000, Q=0.1)
        plt.figure()
        f_raw, p_raw = welch(averaged_spikes_filtered[i,start_time:end_time], fs=10000 )
        plt.semilogy(f_raw, p_raw, label='after filtering')
        plt.xlabel('Frequency (Hz)')
        plt.ylabel('spectral energy density')
        print(f_raw[np.argmax(p_raw)])
        plt.figure()
        plt.plot(t_spike, averaged_spikes_filtered[i,start_time:end_time])
        plt.title(['channel'+ str(i) +' averaged after filtering'])
        plt.show()
### check if common mode can be used for fluctuation denoising

#%%

def cm_noise(averaged_spikes, plot=False,average_xlim = [0, 0.2]):
    averaged_spikes_denoise = copy.deepcopy(averaged_spikes)
    averaged_spikes_denoise[:,1100:1400] = np.squeeze(common_reference_matrix(averaged_spikes[:,1100:1400], ref_type = 'average'))
    if plot:
        t_spike = arange(10000) / 10000
        for i in range(averaged_spikes.shape[0]):
            plt.figure()
            plt.plot(t_spike, averaged_spikes[i,:])
            plt.title(['channel'+ str(i) +' averaged'])
            plt.xlim(average_xlim)
            plt.show()
            plt.figure()
            plt.plot(t_spike, averaged_spikes_denoise[i,:])
            plt.title(['channel'+ str(i) +' averaged after deleting common mode'])
            plt.xlim(average_xlim)
            plt.show()

from scipy.signal import savgol_filter
def smooth_savgol(averaged_spikes,sampling_rate_new,ch_num = None,width = 15, order = 2):
    start_time = 0
    end_time = averaged_spikes.shape[1]
    t_spike = arange(start_time,end_time) / sampling_rate_new
    averaged_spikes_filtered = np.zeros_like(averaged_spikes)
    for i in range(averaged_spikes_filtered.shape[0]):
        averaged_spikes_filtered[i, :] = savgol_filter(averaged_spikes[i, :], width, order)
    if ch_num:
        for i in range(ch_num):
            plt.figure()
            plt.plot(t_spike, averaged_spikes[i,start_time:end_time])
            plt.title(['channel'+ str(i) +' averaged before smoothing'])
            plt.show()
            plt.figure()
            plt.plot(t_spike, averaged_spikes_filtered[i,start_time:end_time])
            plt.title(['channel'+ str(i) +' averaged after smoothing'])
            plt.show()
    return averaged_spikes_filtered


def smooth(x, window_len=11, window='hanning'):
    """smooth the data using a window with requested size.

    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.

    input:
        x: the input signal
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal

    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)

    see also:

    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter

    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """
    #
    # if x.ndim != 1:
    #     raise ValueError, "smooth only accepts 1 dimension arrays."
    #
    # if x.size < window_len:
    #     raise ValueError, "Input vector needs to be bigger than window size."

    if window_len < 3:
        return x

    # if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
    #     raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"

    #s = np.r_[x[window_len - 1:0:-1], x, x[-2:-window_len - 1:-1]]
    # print(len(s))
    if window == 'flat':  # moving average
        w = np.ones(window_len, 'd')
    else:
        w = eval('np.' + window + '(window_len)')

    y = np.convolve(w / w.sum(), x, mode='same')
    return y


def window_smooth(averaged_spikes, sampling_rate_new,ch_num=None, window_len=11, window='hanning'):
    start_time = 0
    end_time = averaged_spikes.shape[1]
    t_spike = arange(start_time,end_time) / sampling_rate_new
    averaged_spikes_filtered = np.zeros_like(averaged_spikes)
    for i in range(averaged_spikes_filtered.shape[0]):
        averaged_spikes_filtered[i, :] = smooth(averaged_spikes[i, :], window_len=window_len, window=window)
    if ch_num:
        for i in range(ch_num):
            plt.figure()
            plt.plot(t_spike, averaged_spikes[i,start_time:end_time])
            plt.title(['channel'+ str(i) +' averaged before '+window+' smoothing'])
            plt.show()
            plt.figure()
            plt.plot(t_spike, averaged_spikes_filtered[i,start_time:end_time])
            plt.title(['channel'+ str(i) +' averaged after '+window+' smoothing'])
            plt.show()
    return averaged_spikes_filtered
from sklearn.decomposition import PCA
def pca_mapping(averaged_spikes_features,n_comp=50,vr_plot=False):
    pca = PCA(n_components=n_comp)
    averaged_spikes_features = np.nan_to_num(averaged_spikes_features)
    pca.fit(averaged_spikes_features)
    variance = pca.explained_variance_
    averaged_spikes_features = pca.transform(averaged_spikes_features)
    if vr_plot:
        plt.figure()
        plt.plot(range(n_comp), variance / sum(variance))
        plt.xlabel('PCs')
        plt.ylabel('variance explained')
    return averaged_spikes_features,variance/sum(variance)

import pickle
def save_obj(obj, name ):
    with open(name , 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)
def load_obj(name ):
    with open(name, 'rb') as f:
        return pickle.load(f)


from skimage import io
import numpy as np
from skimage.morphology import reconstruction
import cv2
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import scipy.linalg
from tifffile import imsave

def cell_recorded(position, device_file, device_raw_file, cell_raw_file,
                  cell_meta, cell_mask_file, gene_type_cluster,cardiomypcyte_cluster,
                  radius,search_method='surface',saving_file=None,plot=False,
                  threshold=50, blockSize=201,near=True):
    # io for reading device image
    device = io.imread(device_file)
    device_raw = io.imread(device_raw_file)
    uplimit = device_raw.shape[0]
    # thresholding to eliminate the non-sensor noise
    device_raw_mp = np.amax(device_raw, axis=0)
    device_raw_mp_filtered = np.interp(device_raw_mp, (device_raw_mp.min(), device_raw_mp.max()), (0, 255))
    device_raw_mp_filtered[device_raw_mp_filtered > threshold] = 0
    # fill holes
    seed = np.copy(device_raw_mp_filtered)
    seed[1:-1, 1:-1] = device_raw_mp_filtered.max()
    mask = device_raw_mp_filtered

    device_raw_mp_filtered = reconstruction(seed, mask, method='erosion')
    device_raw_mp_filtered[device_raw_mp_filtered > 0] = 255

    # Gaussian thresholding to eliminate the non-sensor lines
    gray = device_raw_mp_filtered.astype('uint8')

    # Apply adextract_all_spikesaptive threshold with gaussian size 51x51
    thresh_gray = cv2.adaptiveThreshold(gray, 255, adaptiveMethod=cv2.ADAPTIVE_THRESH_GAUSSIAN_C,
                                        thresholdType=cv2.THRESH_BINARY, blockSize=blockSize, C=0)

    # Find connected components (clusters)
    nlabel, labels, stats, centroids = cv2.connectedComponentsWithStats(thresh_gray, connectivity=8)
    # Find second largest cluster (the cluster is the background):
    sensor_sizes = stats[1:, cv2.CC_STAT_AREA]
    sensor_sizes = sensor_sizes[sensor_sizes > 300]
    var_small = 10000
    for max_size in sensor_sizes:
        max_size_idx = np.where(stats[:, cv2.CC_STAT_AREA] == max_size)[0][0]

        mask_candidate = np.zeros_like(thresh_gray)
        # Draw the cluster on mask
        mask_candidate[labels == max_size_idx] = 255

        # Use "open" morphological operation for removing some artifacts
        mask_candidate = cv2.morphologyEx(mask_candidate, cv2.MORPH_OPEN,
                                          cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (5, 5)))
        #plt.figure()
        if len(sensor_sizes) == 1:
            mask = mask_candidate
        else:
            ret, thresh = cv2.threshold(mask_candidate, 50, 255, cv2.THRESH_BINARY)
            # Find contours:
            contours ,_ = cv2.findContours(thresh, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
            # find centroid
            M = cv2.moments(contours[0])
            mask_cx = round(M['m10'] / M['m00'])
            mask_cy = round(M['m01'] / M['m00'])
            #plt.scatter(mask_cx, mask_cy)
            contours = np.array(contours[0])
            contours = contours.squeeze()
            #plt.scatter(contours[:, 0], contours[:, 1], c='pink', s=0.5)
            var = np.var(np.sqrt(np.sum(np.square(contours - [mask_cx, mask_cy]), axis=1)))
            if var < var_small and len(contours)>50:
                var_small = var
                mask = mask_candidate
            #print(var)
        #plt.imshow(mask, cmap='gray')

    # Use "open" morphological operation for removing some artifacts
    mask = cv2.morphologyEx(mask, cv2.MORPH_OPEN, cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (5, 5)))

    device_position = np.interp(device, (device.min(), device.max()), (0, 255))
    device_position[:, mask == 0] = 0
    device_position[device_position < 0.5 * device_position.max()] = 0

    # use a 2D surface fitting to find the sensor surface
    z, x, y = device_position.nonzero()
    data = np.c_[x, y, z]

    # best-fit linear plane (1st-order)
    A = np.c_[data[:, 0], data[:, 1], np.ones(data.shape[0])]
    C, _, _, _ = scipy.linalg.lstsq(A, data[:, 2])  # coefficients

    # read the cell raw data
    cell_raw = io.imread(cell_raw_file)
    cell_raw_mp = np.amax(cell_raw, axis=0)

    cell_centroid_position = cell_meta[cell_meta['position'] == position][['x', 'y']]
    # read the segmantation mask
    cell_mask = io.imread(cell_mask_file)
    cell_mask_mp = np.amax(cell_mask, axis=0)

    cell_like_centroid = np.zeros_like(cell_raw_mp)
    cell_like_centroid[cell_centroid_position['x'], cell_centroid_position['y']] = 255
    cell_like_centroid_sensor = np.logical_and(mask.transpose(), cell_like_centroid)
    cell_centroid_sensor_position_x, cell_centroid_sensor_position_y = np.nonzero(cell_like_centroid_sensor)
    cell_centroid_position = cell_meta[cell_meta['position'] == position][['z', 'x', 'y']]
    cell_recorded_types = []
    if search_method == 'surface':
        X,Y = mask.nonzero()
        Z = C[0] * X + C[1] * Y + C[2]
        sensor_fit_position = np.c_[Z.astype(int), X.astype(int),Y.astype(int)]
        bad_fit = sensor_fit_position[:,0]>=uplimit
        sensor_fit_position = np.delete(sensor_fit_position, bad_fit, 0)
        sensor_2D_fit = np.zeros_like(cell_mask)
        sensor_2D_fit[sensor_fit_position[:,0],sensor_fit_position[:,1],sensor_fit_position[:,2]] = 255
        contact_cell = np.logical_and(sensor_2D_fit, cell_mask)
        cell_recorded = np.unique(cell_mask[contact_cell]) - 1
        cell_recorded =c
        ell_recorded.tolist()
    elif search_method == 'centroid':
        cell_recorded = []
        cell_recorded_centroid = []
        cell_near = []
        cell_near_centroid = []
        current_near = 100
        for x, y in zip(cell_centroid_sensor_position_x, cell_centroid_sensor_position_y):
            cell_candidates = cell_centroid_position.loc[(cell_centroid_position['x'] == x)
                                                         & (cell_centroid_position['y'] == y)]
            z_recorded = C[0] * cell_candidates['x'].to_numpy() + C[1] * cell_candidates['y'].to_numpy() + C[2]
            #print(cell_candidates['x'].to_numpy(), cell_candidates['y'].to_numpy(), cell_candidates['z'].to_numpy(),
                  #z_recorded)

            distance = (z_recorded - cell_candidates['z'].to_numpy()) * 0.56
            #print(distance)
            if 10 > distance > 0:
                cell_recorded.append(cell_meta.iloc[cell_candidates.index.tolist()[0]].iat[0])
                cell_recorded_centroid.append(cell_candidates)
            if distance < current_near:
                current_near = distance
                cell_near = cell_meta.iloc[cell_candidates.index.tolist()[0]].iat[0]
                cell_near_centroid = cell_candidates
        if (not cell_recorded) and cell_near:
            cell_recorded.append(cell_near)
            cell_recorded_centroid.append(cell_near_centroid)
    elif search_method == 'centroid_nearest':
        celltypes = {'idx': [], 'celltype': []}
        celltype_=[]
        # cell_recorded = []
        # cell_recorded_centroid = []
        cell_near = []
        # cell_near_centroid = []
        current_near = 100
        for x, y in zip(cell_centroid_sensor_position_x, cell_centroid_sensor_position_y):
            cell_candidates = cell_centroid_position.loc[(cell_centroid_position['x'] == x)
                                                         & (cell_centroid_position['y'] == y)]
            z_recorded = C[0] * cell_candidates['x'].to_numpy() + C[1] * cell_candidates['y'].to_numpy() + C[2]
            #print(cell_candidates['x'].to_numpy(), cell_candidates['y'].to_numpy(), cell_candidates['z'].to_numpy(),
                  #z_recorded)

            distance = (z_recorded - cell_candidates['z'].to_numpy()) * 0.56
            #print(distance)
            if 0 <= distance < current_near:
                celltype_ = gene_type_cluster['type'][gene_type_cluster['idx'] == cell_candidates.index.to_numpy()[0]].to_numpy()
                if celltype_.size==0:
                    continue
                else:
                    celltype_ = celltype_[0]
                if celltype_  in cardiomypcyte_cluster:
                    current_near = distance
                    cell_near = cell_meta.iloc[cell_candidates.index.tolist()[0]].iat[0]+1
                    cell_type_near = celltype_
                    # cell_near_centroid = cell_candidates
                else:
                    celltype_=[]
        if cell_near:
            celltypes['idx'].append(cell_near-1)
            celltypes['celltype'].append(cell_type_near)
        if celltypes['celltype']:
            celltypes = pd.DataFrame.from_dict(celltypes)
            cell_recorded = celltypes['idx'].tolist()
            cell_recorded_types = celltypes['celltype'].tolist()
        else:
            cell_recorded = []
            cell_recorded_types = []
        #if (not cell_recorded) and cell_near:
        # cell_recorded.append(cell_near)
        # cell_recorded_centroid.append(cell_near_centroid)
    elif search_method == 'centroid_nearest_revserse':
        celltypes = {'idx': [], 'celltype': []}
        celltype_=[]
        # cell_recorded = []
        # cell_recorded_centroid = []
        cell_near = []
        # cell_near_centroid = []
        current_near = 100
        for x, y in zip(cell_centroid_sensor_position_x, cell_centroid_sensor_position_y):
            cell_candidates = cell_centroid_position.loc[(cell_centroid_position['x'] == x)
                                                         & (cell_centroid_position['y'] == y)]
            z_recorded = C[0] * cell_candidates['x'].to_numpy() + C[1] * cell_candidates['y'].to_numpy() + C[2]
            #print(cell_candidates['x'].to_numpy(), cell_candidates['y'].to_numpy(), cell_candidates['z'].to_numpy(),
                  #z_recorded)

            distance = (cell_candidates['z'].to_numpy() - z_recorded) * 0.56
            #print(distance)
            if 0 <= distance < current_near:
                celltype_ = gene_type_cluster['type'][gene_type_cluster['idx'] == cell_candidates.index.to_numpy()[0]].to_numpy()
                if celltype_.size==0:
                    continue
                else:
                    celltype_ = celltype_[0]
                if celltype_  in cardiomypcyte_cluster:
                    current_near = distance
                    cell_near = cell_meta.iloc[cell_candidates.index.tolist()[0]].iat[0]+1
                    cell_type_near = celltype_
                    # cell_near_centroid = cell_candidates
                else:
                    celltype_=[]
        if cell_near:
            celltypes['idx'].append(cell_near-1)
            celltypes['celltype'].append(cell_type_near)
        if celltypes['celltype']:
            celltypes = pd.DataFrame.from_dict(celltypes)
            cell_recorded = celltypes['idx'].tolist()
            cell_recorded_types = celltypes['celltype'].tolist()
        else:
            cell_recorded = []
            cell_recorded_types = []
        #if (not cell_recorded) and cell_near:
        # cell_recorded.append(cell_near)
        # cell_recorded_centroid.append(cell_near_centroid)
        
        
        
    elif search_method == 'centroid_radius':
        celltypes = {'idx':[],'celltype':[]}
        X, Y = mask.nonzero()
        Z = C[0] * X + C[1] * Y + C[2]
        sensor_fit_position = np.c_[Z.astype(int), X.astype(int), Y.astype(int)]
        bad_fit = np.where(sensor_fit_position[:, 0] >= uplimit)
        sensor_fit_position = np.delete(sensor_fit_position, bad_fit, 0)
        sensor_2D_fit = np.zeros_like(cell_mask)
        sensor_centroid = np.median(sensor_fit_position,axis=0).astype(int)
        sensor_2D_fit[sensor_centroid[0]-radius[0]:sensor_centroid[0]+radius[0],
                        sensor_centroid[1]-radius[1]:sensor_centroid[1]+radius[1],
                        sensor_centroid[2]-radius[2]:sensor_centroid[2]+radius[2]] = 255
        contact_cell = np.logical_and(sensor_2D_fit, cell_mask)
        cell_candidates = np.unique(cell_mask[contact_cell]) - 1
        for sub_cell_id in cell_candidates:
            position_idx = np.array(cell_meta['position'] == position)
            cell_id_idx = np.array(cell_meta['Unnamed: 0'] == sub_cell_id)
            cell_idx = np.logical_and(position_idx, cell_id_idx)
            idx = cell_meta.index[cell_idx].tolist()

            celltype_ = gene_type_cluster['type'][gene_type_cluster['idx'] == idx[0]].to_numpy()
            if celltype_.size==0:
                continue
            else:
                celltype_ = celltype_[0]
            if celltype_  in cardiomypcyte_cluster:
                celltypes['idx'].append(sub_cell_id)
                celltypes['celltype'].append(celltype_)
        #contact_cell = np.logical_and(sensor_2D_fit, cell_mask)
        if celltypes['celltype']:
            celltypes = pd.DataFrame.from_dict(celltypes)
            cell_recorded = celltypes['idx'].tolist()
            cell_recorded_types = celltypes['celltype'].tolist()
        else:
            cell_recorded = []
            cell_recorded_types = []

    elif search_method == 'surface_radius':
        celltypes = {'idx':[],'celltype':[]}
        X, Y = mask.nonzero()
        Z = C[0] * X + C[1] * Y + C[2]
        sensor_fit_position = np.c_[Z.astype(int), X.astype(int), Y.astype(int)]
        bad_fit = np.where(sensor_fit_position[:, 0] >= uplimit)
        sensor_fit_position = np.delete(sensor_fit_position, bad_fit, 0)
        sensor_2D_fit = np.zeros_like(cell_mask)
        sensor_min = np.amin(sensor_fit_position,axis=0).astype(int)
        sensor_max = np.amax(sensor_fit_position,axis=0).astype(int)
        
        sensor_2D_fit[sensor_min[0]-radius[0]:sensor_max[0]+radius[0],
                        sensor_min[1]-radius[1]:sensor_max[1]+radius[1],
                        sensor_min[2]-radius[2]:sensor_max[2]+radius[2]] = 255
        contact_cell = np.logical_and(sensor_2D_fit, cell_mask)
        cell_candidates = np.unique(cell_mask[contact_cell]) - 1
        for sub_cell_id in cell_candidates:
            position_idx = np.array(cell_meta['position'] == position)
            cell_id_idx = np.array(cell_meta['Unnamed: 0'] == sub_cell_id)
            cell_idx = np.logical_and(position_idx, cell_id_idx)
            idx = cell_meta.index[cell_idx].tolist()

            celltype_ = gene_type_cluster['type'][gene_type_cluster['idx'] == idx[0]].to_numpy()
            if celltype_.size==0:
                continue
            else:
                celltype_ = celltype_[0]
            if celltype_  in cardiomypcyte_cluster:
                celltypes['idx'].append(sub_cell_id)
                celltypes['celltype'].append(celltype_)
        #contact_cell = np.logical_and(sensor_2D_fit, cell_mask)
        if celltypes['celltype']:
            celltypes = pd.DataFrame.from_dict(celltypes)
            cell_recorded = celltypes['idx'].tolist()
            cell_recorded_types = celltypes['celltype'].tolist()
        else:
            cell_recorded = []
            cell_recorded_types = []

    elif search_method == 'surface_nearest':
        cell_recorded = []
        cell_recorded_centroid = []
        cell_near = []
        cell_near_centroid = []
        current_near = 100
        X, Y = mask.nonzero()
        Z = C[0] * X + C[1] * Y + C[2]
        sensor_fit_position = np.c_[Z.astype(int), X.astype(int),Y.astype(int)]
        bad_fit = sensor_fit_position[:,0]>=uplimit
        sensor_fit_position = np.delete(sensor_fit_position, bad_fit, 0)
        sensor_2D_fit = np.zeros_like(cell_mask)
        sensor_2D_fit[sensor_fit_position[:,0],sensor_fit_position[:,1],sensor_fit_position[:,2]] = 255
        contact_cell = np.logical_and(sensor_2D_fit, cell_mask)
        cell_recorded_candidates = np.unique(cell_mask[contact_cell]) - 1
        for candidate in cell_recorded_candidates:
            z_recorded = C[0] * cell_centroid_position['x'].iloc[candidate] + C[1] * cell_centroid_position['y'].iloc[candidate]+ C[2]
            distance = (z_recorded - cell_centroid_position['z'].iloc[candidate]) * 0.56
            # print(distance)
            if 0 <= distance < current_near:
                current_near = distance
                cell_near = candidate
                cell_near_centroid = cell_centroid_position.iloc[candidate]
        cell_recorded.append(cell_near)
        cell_recorded_centroid.append(cell_near_centroid)
    elif search_method == 'surface_contact':
        cell_recorded = []
        celltypes = {'idx':[],'celltype':[]}
        X, Y = mask.nonzero()
        Z = C[0] * X + C[1] * Y + C[2]
        sensor_fit_position = np.c_[Z.astype(int), X.astype(int), Y.astype(int)]
        bad_fit = np.where(sensor_fit_position[:, 0] >= uplimit)
        sensor_fit_position = np.delete(sensor_fit_position, bad_fit, 0)
        sensor_2D_fit = np.zeros_like(cell_mask)
        sensor_2D_fit[sensor_fit_position[:, 0], sensor_fit_position[:, 1], sensor_fit_position[:, 2]] = 255
        contact_cell = np.logical_and(sensor_2D_fit, cell_mask)
        cell_candidates = np.unique(cell_mask[contact_cell]) - 1
        for sub_cell_id in cell_candidates:
            position_idx = np.array(cell_meta['position'] == position)
            cell_id_idx = np.array(cell_meta['Unnamed: 0'] == sub_cell_id)
            cell_idx = np.logical_and(position_idx, cell_id_idx)
            idx = cell_meta.index[cell_idx].tolist()

            celltype_ = gene_type_cluster['type'][gene_type_cluster['idx'] == idx[0]].to_numpy()
            if celltype_.size==0:
                continue
            else:
                celltype_ = celltype_[0]
            if celltype_  in cardiomypcyte_cluster:
                celltypes['idx'].append(sub_cell_id)
                celltypes['celltype'].append(celltype_)
        #contact_cell = np.logical_and(sensor_2D_fit, cell_mask)
        if celltypes['celltype']:
            celltypes = pd.DataFrame.from_dict(celltypes)
            cell_recorded = celltypes['idx'].tolist()
            cell_recorded_types = celltypes['celltype'].tolist()
        else:
            cell_recorded = []
            cell_recorded_types = []
        # cell_recorded = []
        # cell_near = []
        # X, Y = mask.nonzero()
        # Z = C[0] * X + C[1] * Y + C[2]
        # sensor_fit_position = np.c_[Z.astype(int), X.astype(int),Y.astype(int)]
        # bad_fit = np.where(sensor_fit_position[:,0]>=66)
        # sensor_fit_position = np.delete(sensor_fit_position, bad_fit, 0)
        # sensor_2D_fit = np.zeros_like(cell_mask)
        # sensor_2D_fit[sensor_fit_position[:,0],sensor_fit_position[:,1],sensor_fit_position[:,2]] = 255
        # contact_cell = np.logical_and(sensor_2D_fit, cell_mask)
        # cell_candidates = np.unique(cell_mask[contact_cell]) - 1
        # for sub_cell_id in cell_candidates:
        #     position_idx = np.array(cell_meta['position'] == position)
        #     cell_id_idx = np.array(cell_meta['Unnamed: 0'] == sub_cell_id)
        #     cell_idx = np.logical_and(position_idx, cell_id_idx)
        #     idx = cell_meta.index[cell_idx].tolist()
        #     celltype = gene_type_cluster['type'][gene_type_cluster['idx']==idx[0]].to_numpy()[0]
        #     if celltype not in cardiomypcyte_cluster:
        #         cell_mask[cell_mask==sub_cell_id+1]=0
        # contact_cell = np.logical_and(sensor_2D_fit, cell_mask)
        # if contact_cell.any():
        #     if position == 36:
        #         print('hello')
        #     if near:
        #         cell_near = np.argmax(np.bincount(cell_mask[contact_cell])) - 1# find the value that appears most of the time
        #     else:
        #         cell_near = np.unique(cell_mask[contact_cell]) -1
        #     cell_recorded.append(cell_near)
        #     if not cell_recorded:
        #         cell_near_centroid = []
        #         current_near = 100
        #         for x, y in zip(cell_centroid_sensor_position_x, cell_centroid_sensor_position_y):
        #             cell_candidates = cell_centroid_position.loc[(cell_centroid_position['x'] == x)
        #                                                          & (cell_centroid_position['y'] == y)]
        #             z_recorded = C[0] * cell_candidates['x'].to_numpy() + C[1] * cell_candidates['y'].to_numpy() + C[2]
        #             distance = (z_recorded - cell_candidates['z'].to_numpy()) * 0.56
        #             if 0 <= distance < current_near:
        #                 current_near = distance
        #                 cell_near = cell_meta.iloc[cell_candidates.index.tolist()[0]].iat[0]
        #                 cell_near_centroid = cell_candidates
        #         cell_recorded.append(cell_near)
        # else:
        #     cell_recorded = []
    elif search_method == 'surface_contact_voting':
        cell_recorded = []
        celltypes = {'idx':[],'celltype':[]}
        X, Y = mask.nonzero()
        Z = C[0] * X + C[1] * Y + C[2]
        sensor_fit_position = np.c_[Z.astype(int), X.astype(int), Y.astype(int)]
        bad_fit = np.where(sensor_fit_position[:, 0] >= uplimit)
        sensor_fit_position = np.delete(sensor_fit_position, bad_fit, 0)
        sensor_2D_fit = np.zeros_like(cell_mask)
        sensor_2D_fit[sensor_fit_position[:, 0], sensor_fit_position[:, 1], sensor_fit_position[:, 2]] = 255
        contact_cell = np.logical_and(sensor_2D_fit, cell_mask)
        cell_candidates = np.unique(cell_mask[contact_cell]) - 1
        for sub_cell_id in cell_candidates:
            position_idx = np.array(cell_meta['position'] == position)
            cell_id_idx = np.array(cell_meta['Unnamed: 0'] == sub_cell_id)
            cell_idx = np.logical_and(position_idx, cell_id_idx)
            idx = cell_meta.index[cell_idx].tolist()

            celltype_ = gene_type_cluster['type'][gene_type_cluster['idx'] == idx[0]].to_numpy()
            if celltype_.size==0:
                continue
            else:
                celltype_ = celltype_[0]
            if celltype_  in cardiomypcyte_cluster:
                celltypes['idx'].append(sub_cell_id)
                celltypes['celltype'].append(celltype_)
        #contact_cell = np.logical_and(sensor_2D_fit, cell_mask)
        if celltypes['celltype']:
            celltypes = pd.DataFrame.from_dict(celltypes)
            #cell_major = np.argmax(np.bincount(celltypes['celltype']))
            cell_major = Counter(celltypes['celltype']).most_common()[0][0]
            celltypes = celltypes[celltypes['celltype'] == cell_major]
            cell_recorded = celltypes['idx'].tolist()
            cell_recorded_types = celltypes['celltype'].tolist()
        else:
            cell_recorded = []
            cell_recorded_types = []
                #cell_recorded_centroid = cell_near_centroid
    elif search_method=='all_position':
        cell_recorded = []
        celltypes = {'idx':[],'celltype':[]}
        cell_candidates = np.unique(cell_mask[cell_mask>0]) - 1
        for sub_cell_id in cell_candidates:
            position_idx = np.array(cell_meta['position'] == position)
            cell_id_idx = np.array(cell_meta['Unnamed: 0'] == sub_cell_id)
            cell_idx = np.logical_and(position_idx, cell_id_idx)
            idx = cell_meta.index[cell_idx].tolist()
            celltype_ = gene_type_cluster['type'][gene_type_cluster['idx'] == idx[0]].to_numpy()
            if celltype_.size==0:
                continue
            else:
                celltype_ = celltype_[0]
            if celltype_  in cardiomypcyte_cluster:
                celltypes['idx'].append(sub_cell_id)
                celltypes['celltype'].append(celltype_)
        #contact_cell = np.logical_and(sensor_2D_fit, cell_mask)
        if celltypes['celltype']:
            celltypes = pd.DataFrame.from_dict(celltypes)
            #cell_major = np.argmax(np.bincount(celltypes['celltype']))
            cell_major = Counter(celltypes['celltype']).most_common()[0][0]
            celltypes = celltypes[celltypes['celltype'] == cell_major]
            cell_recorded = celltypes['idx'].tolist()
            cell_recorded_types = celltypes['celltype'].tolist()
    # if position==44:
    #     print('hello')
    if plot:
        # find the boundary of the binary image
        # Apply cv2.threshold() to get a binary image
        ret, thresh = cv2.threshold(mask, 50, 255, cv2.THRESH_BINARY)
        # Find contours:
        contours ,_ = cv2.findContours(thresh, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
        contours = np.array(contours[0])
        contours = contours.squeeze()
        plt.scatter(contours[:, 0], contours[:, 1], c='pink', s=0.5)
        plt.imshow(cell_mask_mp, 'gray')
        if cell_recorded:
            cell_recorded_idx = [cell + 1 for cell in cell_recorded]
            cell_recorded_mask = np.reshape(np.in1d(cell_mask, cell_recorded_idx), cell_mask.shape)
            cell_mask[cell_recorded_mask == False] = 0
            cell_recorded_mask_mp = np.amax(cell_mask, axis=0)
            cell_recorded_masked = np.ma.masked_where(cell_recorded_mask_mp == 0, cell_recorded_mask_mp)
            #cell_recorded_masked[np.nonzero(cell_recorded_masked)]=(int(cell_major[-1]))*125
            plt.imshow(cell_recorded_masked, 'spring')#tab20b
    if saving_file:
        imsave(saving_file, cell_mask)
    return cell_recorded,cell_recorded_types

from sklearn.neighbors import NearestNeighbors
def cell_recorded_dpc(position, device_file, device_raw_file,
                      dot_meta_file, gene_type_cluster, cardiomypcyte_cluster,knn = False,
                   plot=False, threshold=50, blockSize=201,z_extent=10,single=True,xy_kernel=10):

    dot_pos_file = dot_meta_file + '_position' + str(position) + '_spot_meta.csv'
    dot_meta_pos = pd.read_csv(dot_pos_file, header=0)
    # io for reading device image
    device = io.imread(device_file)
    device_raw = io.imread(device_raw_file)
    uplimit = device_raw.shape[0]
    # thresholding to eliminate the non-sensor noise
    device_raw_mp = np.amax(device_raw, axis=0)
    device_raw_mp_filtered = np.interp(device_raw_mp, (device_raw_mp.min(), device_raw_mp.max()), (0, 255))
    device_raw_mp_filtered[device_raw_mp_filtered > threshold] = 0
    # fill holes
    seed = np.copy(device_raw_mp_filtered)
    seed[1:-1, 1:-1] = device_raw_mp_filtered.max()
    mask = device_raw_mp_filtered

    device_raw_mp_filtered = reconstruction(seed, mask, method='erosion')
    device_raw_mp_filtered[device_raw_mp_filtered > 0] = 255

    # Gaussian thresholding to eliminate the non-sensor lines
    gray = device_raw_mp_filtered.astype('uint8')

    # Apply adextract_all_spikesaptive threshold with gaussian size 51x51
    thresh_gray = cv2.adaptiveThreshold(gray, 255, adaptiveMethod=cv2.ADAPTIVE_THRESH_GAUSSIAN_C,
                                        thresholdType=cv2.THRESH_BINARY, blockSize=blockSize, C=0)

    # Find connected components (clusters)
    nlabel, labels, stats, centroids = cv2.connectedComponentsWithStats(thresh_gray, connectivity=8)
    # Find second largest cluster (the cluster is the background):
    sensor_sizes = stats[1:, cv2.CC_STAT_AREA]
    sensor_sizes = sensor_sizes[sensor_sizes > 300]
    var_small = 10000
    for max_size in sensor_sizes:
        max_size_idx = np.where(stats[:, cv2.CC_STAT_AREA] == max_size)[0][0]

        mask_candidate = np.zeros_like(thresh_gray)
        # Draw the cluster on mask
        mask_candidate[labels == max_size_idx] = 255

        # Use "open" morphological operation for removing some artifacts
        mask_candidate = cv2.morphologyEx(mask_candidate, cv2.MORPH_OPEN,
                                          cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (5, 5)))
        # plt.figure()
        if len(sensor_sizes) == 1:
            mask = mask_candidate
        else:
            ret, thresh = cv2.threshold(mask_candidate, 50, 255, cv2.THRESH_BINARY)
            # Find contours:
            contours, _ = cv2.findContours(thresh, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
            # find centroid
            M = cv2.moments(contours[0])
            mask_cx = round(M['m10'] / M['m00'])
            mask_cy = round(M['m01'] / M['m00'])
            # plt.scatter(mask_cx, mask_cy)
            contours = np.array(contours[0])
            contours = contours.squeeze()
            # plt.scatter(contours[:, 0], contours[:, 1], c='pink', s=0.5)
            var = np.var(np.sqrt(np.sum(np.square(contours - [mask_cx, mask_cy]), axis=1)))
            if var < var_small and len(contours) > 50:
                var_small = var
                mask = mask_candidate
            # print(var)
        # plt.imshow(mask, cmap='gray')

    # Use "open" morphological operation for removing some artifacts
    mask = cv2.morphologyEx(mask, cv2.MORPH_OPEN, cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (5, 5)))
    kernel = np.ones((xy_kernel,xy_kernel),np.uint8)
    
    mask = cv2.dilate(mask,kernel,iterations = 1)
    device_position = np.interp(device, (device.min(), device.max()), (0, 255))
    device_position[:, mask == 0] = 0
    device_position[device_position < 0.5 * device_position.max()] = 0

    # use a 2D surface fitting to find the sensor surface
    z, x, y = device_position.nonzero()
    data = np.c_[x, y, z]

    # best-fit linear plane (1st-order)
    A = np.c_[data[:, 0], data[:, 1], np.ones(data.shape[0])]
    C, _, _, _ = scipy.linalg.lstsq(A, data[:, 2])  # coefficients

    # dots_position = dot_meta_pos[['spot_location_3', 'spot_location_1', 'spot_location_2']]
    celltypes = pd.DataFrame()
    X, Y = mask.nonzero()
    
    Z = C[0] * X + C[1] * Y + C[2]
    sensor_fit_position = np.c_[Z.astype(int), X.astype(int), Y.astype(int)]
    
    for z_ex in range(1,z_extent+1):
        sensor_fit_position=np.vstack((sensor_fit_position,np.c_[Z.astype(int)-z_ex, X.astype(int), Y.astype(int)]))
    bad_fit = np.where(sensor_fit_position[:, 0] >= uplimit)
    sensor_fit_position = np.delete(sensor_fit_position, bad_fit, 0)
    sensor_fit_position = pd.DataFrame(sensor_fit_position, columns=['spot_location_3', 'spot_location_1', 'spot_location_2']) 
    remaining_cells = gene_type_cluster['cell_id'][gene_type_cluster['position'] == position].values
    
    dot_meta_pos = dot_meta_pos[dot_meta_pos['spot_cell_id'].isin(remaining_cells)]
    
    contact_cell_all = pd.merge(dot_meta_pos, sensor_fit_position, how='inner', on=['spot_location_3', 'spot_location_1', 'spot_location_2'])
    contact_cell_all = contact_cell_all.groupby(['spot_cell_id', 'spot_image_position', 'spot_day']).size().reset_index(
        name='count')
    

    for _, row in contact_cell_all.iterrows():
        sub_cell_id = row['spot_cell_id']
        celltype_row = gene_type_cluster.loc[(gene_type_cluster['position'] == position)
                              & (gene_type_cluster['cell_id'] == sub_cell_id)]
        celltype_row.insert(0, "count", row['count'], True)
        celltype_ = celltype_row['cell type'].values
        if celltype_.size == 0:
            continue
        if celltype_ in cardiomypcyte_cluster:
            celltypes = celltypes.append(celltype_row)
    if (single==True) & (celltypes.shape[0]!=0):
        celltypes = celltypes.sort_values(by=['count'], ascending=False).iloc[0]
    if knn == True:
        if len(celltypes)==0:
            cur_len = 0
            print(position)
            sensor_fit_position_mean = np.mean(sensor_fit_position.to_numpy(),axis=0)
            print(sensor_fit_position_mean)
            cell_pos_file = dot_meta_file + '_position' + str(position) + '_cell_meta.csv'
            cell_meta_pos = pd.read_csv(cell_pos_file, header=0)
            print(cell_meta_pos.shape)
            cell_meta_pos = cell_meta_pos[cell_meta_pos['cell_id'].isin(remaining_cells)]
            print(cell_meta_pos.shape)
            cell_centroid = cell_meta_pos[['icl_3','icl_1','icl_2']].to_numpy()
            if cell_centroid.shape[0]>=20:
                nn=20
            else:
                nn=cell_centroid.shape[0]
            neigh = NearestNeighbors(n_neighbors=nn)    
            neigh.fit(cell_centroid)
            cell_id_1nn = neigh.kneighbors([sensor_fit_position_mean], nn, return_distance=False)[0]
            for idd in cell_id_1nn:
                print(idd,cell_meta_pos.iloc[idd]['cell_id'],cell_meta_pos.iloc[idd]['count'])
                c_id = cell_meta_pos.iloc[idd]['cell_id']
                celltype_row = gene_type_cluster.loc[(gene_type_cluster['position'] == position)
                      & (gene_type_cluster['cell_id'] == c_id)]
                celltype_ = celltype_row['cell type'].values
                if celltype_.size == 0:
                    continue
                if celltype_ in cardiomypcyte_cluster:
                    celltypes = celltype_row.iloc[0]
                    break
                else:
                    celltypes = pd.DataFrame()
            print(celltype_)   
            print(celltypes)
        
    if len(celltypes)==0:
        print('No CM cell find in position'+str(position))
        if len(X)<1500:
            print('Because sensor was not found')
        elif contact_cell_all.shape[0] != 0:
            if celltype_ not in cardiomypcyte_cluster:
                print('Because only'+ str(celltype_)+' are found')
        else:
            print('Because no cells nearby')
    if plot:
        # find the boundary of the binary image
        # Apply cv2.threshold() to get a binary image
        ret, thresh = cv2.threshold(mask, 50, 255, cv2.THRESH_BINARY)
        # Find contours:
        contours, _ = cv2.findContours(thresh, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
        contours = np.array(contours[0])
        contours = contours.squeeze()

        if celltypes.shape[0]:
            plt.scatter(dot_meta_pos['spot_location_1'].to_numpy(), dot_meta_pos['spot_location_2'].to_numpy(),
                        c=dot_meta_pos['spot_cell_id'].to_numpy())
            # if celltypes.shape:
            plt.scatter(dot_meta_pos[dot_meta_pos['spot_cell_id']==celltypes['cell_id']]['spot_location_1'].to_numpy(),
                        dot_meta_pos[dot_meta_pos['spot_cell_id']==celltypes['cell_id']]['spot_location_2'].to_numpy(), c='r')
        plt.scatter(contours[:, 0], contours[:, 1], c='pink', s=0.5)

    return celltypes
def contact_cell_finding_dpc(data_folder,save_pkl,dot_meta_file,plot,
                             gene_type_cluster,cardiomypcyte_cluster,single,z_extent=30,xy_kernel=10,threshold=50,ETP=None):
    stage_dir_name = os.walk(data_folder)
    celltypes=pd.DataFrame()
    for path,dir,_ in stage_dir_name:
        for sub_dir in dir:
            if 'position' in sub_dir:
                position = int(sub_dir.split('position')[1])
                if int(position)<=64:
                    sub_dir_name = os.listdir(os.path.join(path,sub_dir))
                    for file in sub_dir_name:
                        if file.endswith('5.tif'):
                            device_file = os.path.join(path,sub_dir,file)
                        elif file.endswith('6.tif'):
                            device_raw_file = os.path.join(path,sub_dir,file)
                    cardiomypcyte_cluster = ETP[position]
                    if plot:
                        plt.figure()
                    celltypes_pos = cell_recorded_dpc(position, device_file, device_raw_file,
                                      dot_meta_file, gene_type_cluster, cardiomypcyte_cluster,
                                     plot=plot,single=single,z_extent=z_extent,xy_kernel=xy_kernel,threshold=threshold)
                    celltypes_pos['position']=position
                    celltypes = celltypes.append(celltypes_pos)
                    if plot:
                        plt.title(['position ',str(position)])
                        plt.show()
    if save_pkl:
        save_obj(celltypes,save_pkl)
    return celltypes

def contact_cell_finding(data_folder,cell_mask_folder,meta_file,save_pkl,plot,
                         search_method,gene_type_cluster,cardiomypcyte_cluster,radius):
    cr_cell_recorded = {'position':[],'cell_id':[],'cell_type':[]}
    stage_dir_name = os.walk(data_folder)
    cell_meta = pd.read_csv(meta_file, header=0)
    for path,dir,_ in stage_dir_name:
        for sub_dir in dir:
            if 'position' in sub_dir:
                position = int(sub_dir.split('position')[1])
                if int(position)<=64:
                    sub_dir_name = os.listdir(os.path.join(path,sub_dir))
                    for file in sub_dir_name:
                        if file.endswith('5.tif'):
                            device_file = os.path.join(path,sub_dir,file)
                        elif file.endswith('6.tif'):
                            device_raw_file = os.path.join(path,sub_dir,file)
                        elif file.endswith('4.tif'):
                            cell_raw_file = os.path.join(path,sub_dir, file)
                    cell_mask_file = os.path.join(cell_mask_folder,'position'+str(position)+'.tif')
                    if plot:
                        plt.figure()

                    cell_id,cell_recorded_types = cell_recorded(position, device_file, device_raw_file, cell_raw_file,
                                     cell_meta, cell_mask_file, gene_type_cluster,cardiomypcyte_cluster,
                                    radius, plot=plot, search_method=search_method,
                                    saving_file=None,threshold=50, blockSize=201)
                    cr_cell_recorded['cell_id'].append(cell_id)
                    cr_cell_recorded['cell_type'].append(cell_recorded_types)
                    cr_cell_recorded['position'].append([position]*len(cell_id))
                    if plot:
                        plt.title(['position ',str(position)])
                        plt.show()
    if save_pkl:
        save_obj(cr_cell_recorded,save_pkl)
    return cr_cell_recorded

def average_bin(arr, n):
    end =  n * int(len(arr)/n)
    return np.mean(arr[:end].reshape(-1, n), 1)
def std_bin(arr, n):
    end =  n * int(len(arr)/n)
    return np.std(arr[:end].reshape(-1, n), 1)
def max_bin(arr, n):
    end =  n * int(len(arr)/n)
    return np.amax(arr[:end].reshape(-1, n), 1)

from collections import deque
def downsample(factor,values):
    buffer_ = deque([],maxlen=factor)
    downsampled_values = []
    print()
    thres = np.std(values)
    for i,value in enumerate(values):
        buffer_.appendleft(value)
        if (i-1)%factor==0:
            #Take max value out of buffer
            max_v = max(buffer_)
            min_v = min(buffer_)
            if max_v-min_v>thres:
                if buffer_.index(max_v)>buffer_.index(min_v):
                    downsampled_values.append(max_v)
                    downsampled_values.append(min_v)
                else:
                    downsampled_values.append(min_v)
                    downsampled_values.append(max_v)
            else:
                mean_v = (max_v+min_v)/2
                downsampled_values.append(mean_v)
                downsampled_values.append(mean_v)
        
    return np.array(downsampled_values)
def downsample_abs(factor,values):
    buffer_ = deque([],maxlen=factor)
    downsampled_values = []
    print()
    thres = np.std(values)
    for i,value in enumerate(values):
        buffer_.appendleft(value)
        if (i-1)%factor==0:
            #Take max value out of buffer
            max_v = max(buffer_)
            min_v = min(buffer_)
#             if max_v-min_v>thres:
            med_v = median(buffer_)
            mean_v = sum(buffer_)/len(buffer_)
            if buffer_.index(max_v)>buffer_.index(min_v):
                downsampled_values.append(max_v)
                if mean_v>med_v:
                    downsampled_values.append(mean_v)
                    downsampled_values.append(med_v)
                    
                else:
                    downsampled_values.append(mean_v)
                    downsampled_values.append(med_v)    
                downsampled_values.append(min_v)
            else:
                downsampled_values.append(min_v)
                downsampled_values.append(med_v)
                downsampled_values.append(max_v)
#             else:
                
#                 mean_v = (max_v+min_v)/2
#                 downsampled_values.append(mean_v) 
#                 downsampled_values.append(mean_v) 
#                 downsampled_values.append(mean_v) 
    return np.array(downsampled_values)
from scipy.stats.mstats import zscore
def plot_heatmap(input_data, labels, feature_name,ax=None, use_imshow=False, show_vlines=True, fontsize=10, annotation=None):
    cmap = plt.cm.get_cmap('bwr')

    if ax is None:
        ax = plt.axes()

    clust_sizes = [sum(labels == i) for i in np.unique(labels)]
    data = np.vstack([input_data.iloc[labels == i, :].loc[:, feature_name].values for i in np.unique(labels)]).T

    if use_imshow:
        ax.imshow(np.flipud(zscore(data, axis=1)), vmin=-2.5, vmax=2.5, cmap=cmap, interpolation='none', aspect='auto')
    else:
        ax.pcolor(np.flipud(zscore(data, axis=1)), vmin=-2.5, vmax=2.5, cmap=cmap)
    # plt.imshow(np.flipud(zscore(data,axis=1)),vmin=-2.5, vmax=2.5, cmap=cmap, aspect='auto', interpolation='none')
    ax.set_xlim([0, data.shape[1]])
    ax.set_ylim([0, data.shape[0]])
    if show_vlines:
        for i in np.cumsum(clust_sizes[:-1]):
            ax.axvline(i, color='k', linestyle='-')
    ax.set_yticks(np.arange(data.shape[0]) + 0.5)

    if not annotation:
        y_labels = feature_name[::-1]
    else:
        y_labels = []
        for gene in feature_name[::-1]:
            y_labels.append("%s - %s" % (gene, annotation[gene]))

    ax.set_yticklabels(y_labels, fontsize=fontsize)
    # ax.get_xaxis().set_fontsize(fontsize)
    plt.tick_params(axis='both', which='major', labelsize=fontsize)


# %%

def plot_heatmap_with_labels(data, feature_name, sample_name, labels, cmap=plt.cm.get_cmap('jet'), show_axis=True,
                             show_top_ticks=True, use_labels=None, font_size=15, annotation=None):
    input_data = pd.DataFrame(data,
                              columns=feature_name, index=sample_name)
    g = plt.GridSpec(2, 1, wspace=0.01, hspace=0.01, height_ratios=[0.5, 10])
    cluster_array = np.expand_dims(np.sort(labels), 1).T
    ax = plt.subplot(g[0])
    ax.imshow(cluster_array, aspect='auto', interpolation='none', cmap=cmap)
    if show_top_ticks:
        locations = []
        for i in np.unique(labels):
            locs = np.median(np.argwhere(cluster_array == i)[:, 1].flatten())
            locations.append(locs)
        ax.xaxis.tick_top()
        if use_labels is not None:
            plt.xticks(locations, use_labels)
        else:
            plt.xticks(locations, np.unique(labels))
        ax.get_yaxis().set_visible(False)
    # ax.axis('off')

    ax = plt.subplot(g[1])
    plot_heatmap(input_data, labels,feature_name, ax=ax, use_imshow=False, show_vlines=True, fontsize=font_size,
                 annotation=annotation)
    if not show_axis:
        plt.axis('off')

def spike_cluster_plot(averaged_spikes, labels, label_name, sampling_rate_new, plot_all=False, save_name=None):
    t = np.arange(averaged_spikes.shape[1]) / sampling_rate_new
    padding = int(averaged_spikes.shape[1]/2)
    for i in np.unique(labels):
        idx = np.array(np.where(labels == i))
        print(idx + 1)
        new_spikes = averaged_spikes[idx, :]
        if new_spikes.ndim == 3:
            new_spikes = new_spikes.squeeze()
        if new_spikes.ndim == 1:
            new_spikes = np.expand_dims(new_spikes, axis=0)

        plt.figure()
        #averaged_spike, new_spikes=spike_aline_average(new_spikes,n_pad_left=padding,n_pad_right=padding,detect_sign = 1)
        t = np.arange(new_spikes.shape[1]) / sampling_rate_new
        averaged_spike = np.mean(new_spikes, axis=0)
        plt.axis('off')
        # new_spikes = wd.align_spikes(new_spikes,new_samp=sampling_rate_new,detect_sign=1)
        averaged_spike_plot(new_spikes, averaged_spike, t)
        plt.title(label_name + str(i))
        if plot_all:
            for j in range(new_spikes.shape[0]):
                plt.figure()
                plt.plot(t, new_spikes[j, :])
                plt.title(label_name + str(i))
        if save_name:
            plt.savefig(save_name + label_name + str(i) + '.pdf'
                        , transparent=True)
            
def spike_aline_average_smooth(spikes,n_pad_left=500,n_pad_right=500,detect_sign = 1):
    length = spikes.shape[1]
    numbers = spikes.shape[0]
#     spikes_new = np.zeros([numbers, n_pad_left + n_pad_right])
    spikes_new = []
    # cur_move=0
    sp_var = np.var(spikes, axis=1)
    flag=0
    for i in range(numbers):
        if sp_var[i]>0.5:
            spikes[i,:] = smooth(spikes[i,:], window_len=20, window='hanning')
            
        if detect_sign == -1:
            peak_idx = np.argmin(spikes[i,:int(length/2+length/10)])
        elif detect_sign == 1:
            peak_idx = np.argmax(spikes[i,:int(length/2)])
        elif detect_sign == 'diff':
            peak_idx_abs = np.argmax(np.abs(spikes[i,:]))

            spikes_diff = np.diff(spikes[i,int(peak_idx_abs-length/20):int(peak_idx_abs+length/20)])
            if len(spikes_diff)>0:
                peak_idx = np.argmin(spikes_diff)+int(peak_idx_abs-length/20)
            else:
                peak_idx = peak_idx_abs
                flag=1
        elif detect_sign == 'abs':
            peak_idx = np.argmax(np.abs(spikes[i,:]))

        if peak_idx - n_pad_left > 0 and peak_idx + n_pad_right < length:
            spike = spikes[i,peak_idx - n_pad_left:peak_idx + n_pad_right]
        elif peak_idx - n_pad_left < 0:
            spike = spikes[i,:peak_idx + n_pad_right]
            spike = np.pad(spike, (np.abs(peak_idx - n_pad_left), 0), 'constant')
        elif peak_idx + n_pad_right > length:
            spike = spikes[i,peak_idx - n_pad_left:]
            spike = np.pad(spike, (0, peak_idx + n_pad_right - length), 'constant')
        if flag==1:
            spike=np.zeros(len(spike))
            flag=0
        # if cur_move<np.abs(peak_idx-int(length/2)):
        #     cur_move=np.abs(peak_idx-int(length/2))
#         spikes_new[i, :] = spike
        spikes_new.append(spike)
    min_len = min([len(ti) for ti in spikes_new])
    spikes_new = np.array([ti[:min_len] for ti in spikes_new])
    # spikes_new=spikes_new[:,cur_move:]
#     spikes_new=spikes_new[:,cur_move:-1-cur_move]
    averaged_spike = np.mean(spikes_new, axis=0)
    return averaged_spike, spikes_new
            
def spike_aline_average(spikes,n_pad_left=500,n_pad_right=500,detect_sign = 1):
    length = spikes.shape[1]
    numbers = spikes.shape[0]
#     spikes_new = np.zeros([numbers, n_pad_left + n_pad_right])
    spikes_new = []
    # cur_move=0
    for i in range(numbers):
        if detect_sign == -1:
            peak_idx = np.argmin(spikes[i,:int(length/2+length/10)])
        elif detect_sign == 1:
            peak_idx = np.argmax(spikes[i,:int(length/2)])
        elif detect_sign == 'diff':
            spikes_diff = np.diff(spikes[i,:])
            peak_idx = np.argmin(spikes_diff)
        elif detect_sign == 'abs':
            peak_idx = np.argmax(np.abs(spikes[i,:]))
            
        if peak_idx - n_pad_left > 0 and peak_idx + n_pad_right < length:
            spike = spikes[i,peak_idx - n_pad_left:peak_idx + n_pad_right]
        elif peak_idx - n_pad_left < 0:
            spike = spikes[i,:peak_idx + n_pad_right]
            spike = np.pad(spike, (np.abs(peak_idx - n_pad_left), 0), 'constant')
        elif peak_idx + n_pad_right > length:
            spike = spikes[i,peak_idx - n_pad_left:]
            spike = np.pad(spike, (0, peak_idx + n_pad_right - length), 'constant')

        # if cur_move<np.abs(peak_idx-int(length/2)):
        #     cur_move=np.abs(peak_idx-int(length/2))
#         spikes_new[i, :] = spike
        spikes_new.append(spike)
    min_len = min([len(ti) for ti in spikes_new])
    spikes_new = np.array([ti[:min_len] for ti in spikes_new])
    # spikes_new=spikes_new[:,cur_move:]
#     spikes_new=spikes_new[:,cur_move:-1-cur_move]
    averaged_spike = np.mean(spikes_new, axis=0)
    return averaged_spike, spikes_new
from matplotlib.ticker import MaxNLocator
def spike_cluster_overlay(averaged_spikes, labels, label_name, sampling_rate_new, title=None, plot_all=False,
                          save_name=None):
    t = np.arange(averaged_spikes.shape[1]) / sampling_rate_new
    ax = plt.figure().gca()
    for i in np.unique(labels):
        idx = np.array(np.where(labels == i))
        print(idx + 1)
        new_spikes = averaged_spikes[idx, :]
        if new_spikes.ndim == 3:
            new_spikes = new_spikes.squeeze()
        if new_spikes.ndim == 1:
            new_spikes = np.expand_dims(new_spikes, axis=0)
        averaged_spike = np.mean(new_spikes, axis=0)
        plt.plot(t, averaged_spike, '.-', label=i)
        plt.legend()
        plt.title(title)
        ax.xaxis.set_major_locator(MaxNLocator(
            integer=True))  # new_spikes = wd.align_spikes(new_spikes,new_samp=sampling_rate_new,detect_sign=1)
        if plot_all:
            for j in range(new_spikes.shape[0]):
                plt.figure()
                plt.plot(t, new_spikes[j, :])
                plt.title(label_name + str(i))

        if save_name:
            plt.savefig(save_name + label_name + str(i) + '.pdf'
                        , transparent=True)

import seaborn as sns
from matplotlib import rcParams
import scanpy as sc
def position_std_plot(seq_stage1, name):
    res = pd.DataFrame(columns=seq_stage1.var_names, index=seq_stage1.obs['position'].unique())
    for pos in seq_stage1.obs['position'].unique():
        res.loc[pos] = seq_stage1[seq_stage1.obs['position'].isin([pos]), :].X.mean(0)
    sns.distplot(res.std(), axlabel=name)
    plt.show()
    print(name)
    print(res.std().nlargest(10))
    return res


def gene_data_collection(seq_stage_files, gene_name_files, meta_files, batch_categories, sen_choose = False,
                         cut_off=[60, 40, 30, 10], pre=False, position_filt=False,sen_wi=True,pos_wi = list(range(1, 65))):
    seq_stage1_file = seq_stage_files[0]
    meta_file = meta_files[0]
    seq_stage1 = sc.read_csv(seq_stage1_file)
    gene_name = pd.read_csv(gene_name_files[0], header=None)
    seq_stage1.var_names = gene_name[2]
    cell_meta = pd.read_csv(meta_file, header=0)
    cell_meta['cell_num'] = cell_meta['Unnamed: 0'].groupby(cell_meta['position']).rank(ascending=1, method='dense')

    seq_stage1.obs['position'] = cell_meta['position'].values
    seq_stage1.obs['cell_num'] = cell_meta['cell_num'].values
    seq_stages = []
    # cut_off=[60,40,30,10]
    #pos_wi = list(range(1, 65))
    if sen_choose:
        if sen_wi:
            seq_stage1 = seq_stage1[seq_stage1.obs['position'].isin(pos_wi)]
        else:
            seq_stage1 = seq_stage1[~seq_stage1.obs['position'].isin(pos_wi)]
    if position_filt:
        seq_stage1.obs['position'] = seq_stage1.obs['position'].astype('category')
        sc.pp.calculate_qc_metrics(seq_stage1, percent_top=None, log1p=False, inplace=True)
        f, ax = plt.subplots(figsize=[15, 1.5])
        sc.pl.violin(seq_stage1, keys=['total_counts'], groupby='position', ax=ax)
        plt.show()
        res_mean = pd.DataFrame(columns=['mean'], index=seq_stage1.obs['position'].unique())
        for pos in seq_stage1.obs['position'].unique():
            res_mean.loc[pos] = np.array(seq_stage1[seq_stage1.obs['position'].isin([pos]), :].X.sum(1)).mean()
        sns.distplot(res_mean)
        plt.show()


    if pre:
        sc.pp.filter_cells(seq_stage1, min_genes=10)  # change to count
        sc.pp.filter_cells(seq_stage1, min_counts=cut_off[0])  # change to count
        sc.pp.calculate_qc_metrics(seq_stage1, percent_top=None, log1p=False, inplace=True)
        sc.pl.violin(seq_stage1, ['n_genes_by_counts', 'total_counts'],
                     jitter=4, multi_panel=True)

        sc.pl.scatter(seq_stage1, x='total_counts', y='n_genes_by_counts')
        sc.pp.normalize_total(seq_stage1, target_sum=150)
        sc.pp.log1p(seq_stage1)
        seq_stage1.raw = seq_stage1
        res = position_std_plot(seq_stage1, batch_categories[0] + ' std between positions')
        sc.pp.combat(seq_stage1, key='position')
        res = position_std_plot(seq_stage1, batch_categories[0] + ' std between positions')
        sc.tl.pca(seq_stage1, svd_solver='arpack')
        sc.pl.pca(seq_stage1, color='position')

    for i in range(1, len(batch_categories)):
        seq_stage_file_ = seq_stage_files[i]
        meta_file = meta_files[i]
        seq_stage_ = sc.read_csv(seq_stage_file_)
        gene_name = pd.read_csv(gene_name_files[i], header=None)
        seq_stage_.var_names = gene_name[2]
        cell_meta = pd.read_csv(meta_file, header=0)
        cell_meta['cell_num'] = cell_meta['Unnamed: 0'].groupby(cell_meta['position']).rank(ascending=1, method='dense')
        seq_stage_.obs['position'] = cell_meta['position'].values
        seq_stage_.obs['cell_num'] = cell_meta['cell_num'].values
        if sen_choose:
            if sen_wi:
                seq_stage_ = seq_stage_[seq_stage_.obs['position'].isin(pos_wi)]
            else:
                seq_stage_ = seq_stage_[~seq_stage_.obs['position'].isin(pos_wi)]
        if position_filt:
            seq_stage_.obs['position'] = seq_stage_.obs['position'].astype('category')
            sc.pp.calculate_qc_metrics(seq_stage_, percent_top=None, log1p=False, inplace=True)
            f, ax = plt.subplots(figsize=[15, 1.5])
            sc.pl.violin(seq_stage_, keys=['total_counts'], groupby='position', ax=ax)
            plt.show()
            res_mean = pd.DataFrame(columns=['mean'], index=seq_stage_.obs['position'].unique())
            for pos in seq_stage_.obs['position'].unique():
                res_mean.loc[pos] = np.array(seq_stage_[seq_stage_.obs['position'].isin([pos]), :].X.sum(1)).mean()
            sns.distplot(res_mean)
            plt.show()

        if pre:
            sc.pp.filter_cells(seq_stage_, min_genes=10)  # change to count
            sc.pp.filter_cells(seq_stage_, min_counts=cut_off[i])  # change to count
            sc.pp.calculate_qc_metrics(seq_stage_, percent_top=None, log1p=False, inplace=True)
            sc.pl.violin(seq_stage_, ['n_genes_by_counts', 'total_counts'],
                         jitter=4, multi_panel=True)
            sc.pl.scatter(seq_stage_, x='total_counts', y='n_genes_by_counts')
            sc.pp.normalize_total(seq_stage_, target_sum=150)
            sc.pp.log1p(seq_stage_)
            seq_stage_.raw = seq_stage_
            res = position_std_plot(seq_stage_, batch_categories[i] + ' std between positions')
            sc.pp.combat(seq_stage_, key='position')
            res = position_std_plot(seq_stage_, batch_categories[i] + ' std between positions')
            sc.tl.pca(seq_stage_, svd_solver='arpack')
            sc.pl.pca(seq_stage_, color='position')

        seq_stages.append(seq_stage_)
    stage_all = seq_stage1.concatenate(seq_stages, batch_key='stage',
                                       batch_categories=batch_categories)
    return stage_all

