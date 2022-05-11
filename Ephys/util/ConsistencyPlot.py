#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from sklearn.metrics import confusion_matrix
import seaborn as sn
import matplotlib.pyplot as plt
import pylab
import matplotlib
from scipy.optimize import linear_sum_assignment
import plotly.graph_objects as go
import plotly.express as pex
import numpy as np
import pandas as pd
def color_cm(cmap,NUM_COLORS,):
    color = []
    color_idx = 0
    cm = pylab.get_cmap(cmap)
    for i in range(NUM_COLORS):
        color.append(matplotlib.colors.to_hex(cm(1. * i / NUM_COLORS)))  # color will now be an RGBA tuple
    return color
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
def matrix_scatterplot(M, xticklabels, yticklabels, xlabel='', ylabel='', fig_width=10, fig_height=14, scale_factor=10.0):
    """Plots a matrix with points as in a scatterplot. Area of points proportional to each matrix element. 
    Suitable to show sparse matrices.

    Args:
        M (np.array): a 2D array
        xticklabels: label list
        yticklabels: label list
        fig_width (int): matplotlib figure width
        fig_height (int): matplotlib figure height
        scale_factor (float): scales the points by this value. 
    """
    Mplot = M.copy()*scale_factor
    Mplot = np.flip(Mplot, axis=0)
    yticklabels.reverse()
    x = np.arange(0, M.shape[1], 1)
    y = np.arange(0, M.shape[0], 1)
    xx, yy = np.meshgrid(x, y)
    plt.figure(figsize=(fig_width, fig_height))
#     plt.scatter(np.ravel(xx), np.ravel(yy), s=np.ravel(Mplot), c='dodgerblue')
    cmap='Set1'
    color = color_cm(cmap,M.shape[0]+M.shape[1])
    colors = np.repeat(color[:M.shape[0]],M.shape[1])
    plt.scatter(np.ravel(xx), np.ravel(yy), s=np.ravel(Mplot), c=colors[::-1])
    
    
    ax = plt.gca()
    ax.set_xlim(np.min(x)-0.5, np.max(x)+0.5)
    ax.set_ylim(np.min(y)-0.5, np.max(y)+0.5)
    ax.set_xticks(x)
    ax.set_xticklabels(xticklabels, rotation=90)
    ax.set_yticks(y)
    ax.set_yticklabels(yticklabels, rotation=0)
    ax.xaxis.set_ticks_position('top')
    ax.yaxis.set_ticks_position('left')
    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)

    ax.tick_params(color='None')
    ax.spines["top"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)

    for tick in ax.get_yticklabels():
        #tick.set_fontname("DejaVu Sans Mono")
        tick.set_fontfamily('monospace')
        tick.set_fontsize(12)

    for tick in ax.get_xticklabels():
        tick.set_fontfamily('monospace')
        tick.set_fontsize(12)

    plt.grid(color='gray', linestyle='-', linewidth=0.5, alpha=0.4)
    plt.box(False)
    plt.tight_layout()
    plt.show()
    return
def consistency_plots(t_type,e_type,t_name='T',e_name='E',re_order=True):

    C_e_types = contingency(a=e_type,
                  b=t_type,
                  unique_a=np.unique(e_type),
                  unique_b=np.unique(t_type))
    #Assign labels of clusters based on 'best match' with transcriptomic celltype label
    row_ind,col_ind = linear_sum_assignment(-C_e_types)

    C_ordered = C_e_types[:,col_ind]
    order_y = np.unique(t_type)[col_ind]

    t_labels_matched = t_type.copy()
    e_labels_matched = e_type.copy()
    if re_order:
        for name in order_y:
            ind = t_type == name
            t_labels_matched[ind] = name

    C_TE_consistency = contingency(a=e_labels_matched,
                    b=t_labels_matched,
                    unique_a=np.unique(e_labels_matched),
                    unique_b=np.unique(t_labels_matched))

    ax_xticklabels = [t_name+'-{}{:>4s}'.format(x,'({:d})'.format(np.sum(t_labels_matched==x))) for x in np.unique(t_labels_matched)]
    ax_yticklabels = [e_name+'-{} ({:d})'.format(y,np.sum(e_labels_matched==y)) for y in np.unique(e_labels_matched)]

    matrix_scatterplot(M=C_TE_consistency,scale_factor=20,
                       xticklabels=ax_xticklabels,
                       yticklabels=ax_yticklabels,
                       xlabel=r'Transcriptomic type',
                       ylabel=r'Electrophysiology type',
                       fig_width=8,fig_height=8)

    ax_xticklabels = [t_name+'-{}'.format(x) for x in np.unique(t_labels_matched)]
    ax_yticklabels = [e_name+'-{}'.format(y) for y in np.unique(e_labels_matched)]
    df_et = pd.DataFrame(C_TE_consistency, columns = ax_xticklabels,
                       index= ax_yticklabels)
    plt.figure(figsize = (10,7))
    sn.heatmap(df_et, annot=True)
    plt.title('gene clustering VS ephys cluster ')
    plt.savefig("stage1+2+3_cm.pdf", transparent=True)
    plt.figure()
    sn.heatmap(df_et.div(df_et.sum(axis=1), axis=0), annot=True)
    plt.title('gene clustering VS  ephys cluster ')
    # cm_1 = confusion_matrix(t_type, e_type, labels=['CM 1','CM 2','CM 3'],normalize='pred')
    # df_cm = pd.DataFrame(cm_1, index = ['gene CM '+i for i in "123"],
    #                   columns = ['ephys CM '+i for i in "123"])
    # plt.figure(figsize = (10,7))
    # sn.heatmap(df_cm, annot=True)
    # plt.title('gene clustering VS GMM ephys cluster ')

    # plt.savefig("stage1+2+3_cmn.pdf", transparent=True)    
    df_river = {'e type':[],'t type':[],'et value':[]}
    for t_cell in df_et.columns.to_list():
        for e_cell in df_et.index.to_list():
            df_river['e type'].append(e_cell)
            df_river['t type'].append(t_cell)
            df_river['et value'].append(df_et.loc[e_cell,t_cell])
    df_river = pd.DataFrame(df_river)

    all_nodes = df_river['e type'].unique().tolist() + df_river['t type'].unique().tolist()
    source_indices = [all_nodes.index(e_type) for e_type in df_river['e type']]
    target_indices = [all_nodes.index(t_type) for t_type in df_river['t type']]

    NUM_COLORS = len(all_nodes)
    cmap='Set1'
    colors = color_cm(cmap,NUM_COLORS)


    # colors = pex.colors.qualitative.D3

    # node_colors_mappings = dict([(node,np.random.choice(colors)) for node in all_nodes])
    node_colors_mappings = dict([(node,colors[idx]) for idx,node in enumerate(all_nodes)])
    node_colors = [node_colors_mappings[node] for node in all_nodes]
    edge_colors = [node_colors_mappings[node] for node in df_river['e type']]

    fig = go.Figure(data=[go.Sankey(
        node = dict(
          pad = 20,
          thickness = 20,
          line = dict(color = "black", width = 1.0),
          label =  all_nodes,
          color =  node_colors,
        ),

        link = dict(
          source =  source_indices,
          target =  target_indices,
          value =  df_river['et value'],
          color = edge_colors,
    ))])

    fig.update_layout(title_text="e type to t type",
                      height=600,
                      font_size=10)
    fig.show()
    
    

