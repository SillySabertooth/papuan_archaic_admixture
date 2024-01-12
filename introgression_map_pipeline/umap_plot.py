import math 
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd 
import seaborn as sns
from sklearn.preprocessing import StandardScaler
#import umap
import umap
import sys

sns.set(rc={'figure.figsize':(15,15)})

#python3 ../arch_pipeline_snakemake/umap_plot.py res_concat/seq_sim.collapsed_annotated.tsv 
#module load python/3.9.12
#source ~/virr_envs/umap_env/bin/activate
#module load py-mpmath/1.2.1 py-numpy/1.22.4 py-pandas/1.4.2 py-seaborn/0.11.2 py-scikit-learn/1.1.1 py-matplotlib/3.5.2 

file_=sys.argv[1]
file_name=file_[:-7]

dt = pd.read_csv(file_, compression='gzip', sep='\t') #read data
dt.replace(np.inf,np.nan, inplace=True)

# default!

dtsbse=dt[dt.ILS != "ILS_FAIL"].reset_index(drop=True)

#dtsbse
X_tsne = dtsbse.iloc[:,[3,4,5,6]]
X_tsne = X_tsne.dropna()


scaler = StandardScaler()
scaler.fit(X_tsne)

X_tsne = scaler.transform(X_tsne)

cols_interest=dtsbse[~dtsbse.iloc[:,[3,4,5,6]].isna().any(axis=1)].reset_index(drop=True)

## umap

reducer = umap.UMAP()
embedding = reducer.fit_transform(X_tsne) 
umapsdf = pd.DataFrame({'umap_dim1':embedding.T[0], 'umap_dim2':embedding.T[1]})

tsnedf_plot = cols_interest.join(umapsdf)


sns.scatterplot(data=tsnedf_plot, x="umap_dim1", y="umap_dim2", hue="closer_arch", alpha = 0.6).set(
    title='umap.ILS_PASS')
plt.savefig(file_name+".full_plot.umap.png", bbox_inches='tight',dpi=300)
plt.clf()

sbset_data = tsnedf_plot[tsnedf_plot.closer_arch.isin(['AltaiNeandertal_closer','Denisova_closer','Vindija33.19_closer','Chagyrskaya-Phalanx_closer','Neand_Denisova','Neand'])]
sns.scatterplot(data=sbset_data, x="umap_dim1", y="umap_dim2", hue="closer_arch", alpha = 0.6).set(
    title='umap.only_introgressed_plotted.ILS_PASS')
plt.savefig(file_name+".introgres_plot.umap.png", bbox_inches='tight',dpi=300)
plt.clf()


# try to plot ILS_FAIL

X_tsne = dt.iloc[:,[3,4,5,6]]
X_tsne = X_tsne.dropna()


scaler = StandardScaler()
scaler.fit(X_tsne)

X_tsne = scaler.transform(X_tsne)

cols_interest=dt[~dt.iloc[:,[3,4,5,6]].isna().any(axis=1)].reset_index(drop=True)

## umap

reducer = umap.UMAP()
embedding = reducer.fit_transform(X_tsne) 
umapsdf = pd.DataFrame({'umap_dim1':embedding.T[0], 'umap_dim2':embedding.T[1]})

tsnedf_plot = cols_interest.join(umapsdf)


sns.scatterplot(data=tsnedf_plot, x="umap_dim1", y="umap_dim2", hue="closer_arch", alpha = 0.6).set(
    title='umap.ILS_ALL')
plt.savefig(file_name+".full_plot.ILS_all.umap.png", bbox_inches='tight',dpi=300)
plt.clf()
