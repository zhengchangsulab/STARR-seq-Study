#!/usr/bin/env python
# coding: utf-8

# In[6]:


#!/usr/bin/python
from __future__ import print_function
import numpy as np
import sys
import matplotlib
#get_ipython().magic(u'matplotlib inline')
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.ticker import StrMethodFormatter, ScalarFormatter
plt.rcParams['svg.fonttype'] = 'none'


cell_name = sys.argv[1]
set_legend = sys.argv[2]
score = sys.argv[3]


#min_value = -8
#max_value = 6

min_value = -4
max_value = 4


def setting_axes(ax, score, set_flag):
    for axis in ['top', 'bottom', 'left', 'right']:
        ax.spines[axis].set_linewidth(4)
        
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    #ax.set_xlim(-9, 9)
    #ax.set_ylim(0, 0.75)
    ax.set_ylim(0, 1.2)
    
    #ax.set_xticks(np.arange(-9.0, 9.0, 2))
    #ax.set_yticks([0, 0.1, 0.2, 0.3])
    #ax.xaxis.set_tick_params(labelsize=30)
    #ax.yaxis.set_tick_params(labelsize=30)
    #ax.set_xticklabels(ax.get_xticks(), fontsize=30, rotation=0)
    #ax.set_yticklabels(ax.get_yticks(), fontsize=30, rotation=40)
    #ax.set_yticklabels([0., 0.05, 0.1,  0.15, 0.2,  0.25, 0.3], fontsize=30, rotation=40)

    #ax.set_yticklabels([str(round(y, 1)) for y in ax.get_yticks()], fontsize=30, rotation=0)
    #ax.set_yticklabels([0., 0.1, 0.2, 0.3], fontsize=30, rotation=0)
    #ax.set_yticklabels(ax.get_yticks(), fontsize=30, rotation=0)
    
    #ax.yaxis.set_major_formatter(StrMethodFormatter('{x:,.2f}'))
    #ax.yaxis.set_major_formatter(StrMethodFormatter('{x:,.4f}'))
    
    xfmt = ScalarFormatter(useMathText=True)
    xfmt.set_powerlimits((-1,2))
    #ax.yaxis.set_major_formatter(xfmt)
    ax.xaxis.set_major_formatter(xfmt)
    
    ax.yaxis.get_offset_text().set_fontsize(30)
    ax.xaxis.get_offset_text().set_fontsize(30)
    
    if score == "gerp":
        score = "GERP"
    else:
        score = "phyloP"

    ax.set_xlabel(score+" "+"score", fontsize=35)
    ax.set_ylabel("Density", fontsize=35)

    ax.tick_params(axis='both', length=6, width=3, labelsize=30)
    
    ax.locator_params(axis='x', nbins=10)
    ax.locator_params(axis='y', nbins=10)
    ax.xaxis.set_label_coords(0.5, -0.15)
    

    ax.legend(fontsize=5, loc='right', bbox_to_anchor=(1.0, 0.9), frameon=False, borderaxespad=0)
    
    if set_flag == "False":
        ax.legend().set_visible(False)


fig, ax = plt.subplots(1,1, figsize=(8,8), dpi=300, facecolor='white', sharex=False, sharey=False)

color = {'a':'g', 'b':'r', 'c':'b', 'd':'c', 'e':'m', 'f':'k'}
label = {'a':'A', 'b':'B', 'c':'C', 'd':'D', 'e':'E', 'f':'F'}
linetypes = {'a':'solid', 'b':'solid', 'c':'solid', 'd':'solid', 'e':'dashed', 'f':'dashed'}

for i in ['a', 'b', 'c', 'd', 'e', 'f']:
    region_name = "{}.region.{}.remove_nan.{}.npy".format(cell_name, i, score)
    region_conservation = np.load(region_name)
    ax = sns.kdeplot(region_conservation, bw=0.3, clip=(min_value, max_value), color=color[i], linewidth=2,linestyle=linetypes[i], label=label[i], ax=ax)


ax.axvline(x=-1, linestyle='dashed', color='r', linewidth=2)
ax.axvline(x=1, linestyle='dashed', color='r', linewidth=2)



setting_axes(ax, score, set_flag=set_legend)

plt.tight_layout()

plt.savefig("Compare_{}_abcedf_{}.{}.drop.nan.density.v2.png".format(cell_name, score, set_legend), dpi=300)
plt.savefig("Compare_{}_abcedf_{}.{}.drop.nan.density.v2.svg".format(cell_name, score, set_legend), dpi=300)

