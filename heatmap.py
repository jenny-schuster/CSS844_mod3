import os
import seaborn as sns
import pandas as pd
import numpy as np
from scipy.stats import zscore
from scipy.cluster.hierarchy import linkage, leaves_list, dendrogram

import matplotlib.pyplot as plt

#Insert path to module/trait correlation matrix
data = pd.read_table("/your/file/here", sep = '\t', header=0, index_col=0)


x_labels = data.columns.values.tolist()

y_labels = list(data.index.values)
print(data)
print(y_labels)
print(x_labels)
for value in range(0, len(y_labels)):
	y_labels[value]=y_labels[value][2:].capitalize()


plt.figure(figsize=(3,6))
hm = sns.heatmap(data, cmap = 'PiYG', xticklabels = x_labels, yticklabels = y_labels, center = 0, clip_on=False, vmin=-0.5, vmax=0.5)
#hm = sns.heatmap(data, cmap = 'PiYG', xticklabels = x_labels, yticklabels = y_labels, center = 0, linewidth=0.25, linecolor="black", cbar_kws={"shrink": 0.75}, clip_on=False, vmin=-0.5, vmax=0.5)
hm.set_xlabel("Trait", fontsize = 'medium', fontweight = 'bold')
hm.set_ylabel("Modules", fontsize = 'medium', fontweight = 'bold')
hm.set_title("Module Trait Correlations", fontsize = 'large', fontweight = 'bold')
plt.yticks(fontsize = 'xx-small')
plt.xticks(fontsize = 'small')

plt.tight_layout()
fig = hm.get_figure()
plt.show()
