import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns
import scipy
import mpmath

from sklearn.svm import SVC
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import MinMaxScaler
from sklearn.preprocessing import RobustScaler
from sklearn.manifold import TSNE
from sklearn.feature_selection import VarianceThreshold

from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.cluster.hierarchy import fcluster

from quantiprot.utils.io import load_fasta_file
from quantiprot.utils.feature import Feature, FeatureSet
from quantiprot.utils.sequence import subset, columns, SequenceSet
from quantiprot.utils.sequence import compact

from quantiprot.utils.mapping import simplify
from quantiprot.metrics.aaindex import get_aaindex_file
from quantiprot.metrics.aaindex import get_aa2charge, get_aa2hydropathy
from quantiprot.metrics.basic import average
from quantiprot.metrics.ngram import pattern_match, pattern_count

from quantiprot.utils.sequence import compact, subset
from quantiprot.metrics.ngram import NgramFeatureSet
from quantiprot.metrics.alphabet import PROTEIN




from Bio import SeqIO
#Load sequence
length_seqs = []
for record in SeqIO.parse("sequence_2.fasta","fasta"):
    length_seqs.append(len(record))
    #print((record))
    

#load the sequence from the file
seq = load_fasta_file("sequence_2.fasta")
SequenceIds = []
SequenceIds2_list = []
for i in SequenceSet.ids(seq):
    SequenceIds.append(i)
for i in SequenceIds:
    SequenceIds2 = i[i.find("[")+1:i.find("]")]
    SequenceIds2_list.append(SequenceIds2)
    

#gather important protein features
polarity = Feature(get_aaindex_file("GRAR740102")).then(average)
hydropathy = Feature(get_aaindex_file("KYTJ820101")).then(average)
iso_point = Feature(get_aaindex_file("ZIMJ680104")).then(average)
pk_COOH = Feature(get_aaindex_file("JOND750102")).then(average)
entropy_form = Feature(get_aaindex_file("HUTJ700103")).then(average)
melting_point = Feature(get_aaindex_file("FASG760102")).then(average)
net_charge = Feature(get_aaindex_file("KLEP840101")).then(average)
glycine = Feature(pattern_count, pattern = 'G')
RGD = Feature(pattern_count, pattern = 'RGD')
GFPGER = Feature(pattern_count, pattern = 'GFPGER')

#Build the feature set
fs=FeatureSet("my set")
fs1 = FeatureSet("test")
fs2 = FeatureSet("glycine")
fs3 = FeatureSet("GFPGER")

#Add the feature to new feature set:
fs.add(polarity, name = "polarity")
fs.add(hydropathy, name = "hydropathy")
fs.add(iso_point, name = "isoelectric point")
fs.add(pk_COOH, name = "pK-COOH")
fs.add(entropy_form, name = "entropy of formation")
fs.add(melting_point, name = "melting point")
fs.add(net_charge, name = "net charge")
fs.add(glycine, name = "glycine")
fs.add(RGD, name = "RGD")
Features = ["polarity", "hydropathy", "isoelectric point", "pK-COOH", "entropy of formation", "melting point", "net charge", "glycine", "RGD"]
Features2 = ["polarity", "hydropathy", "isoelectric point", "pK-COOH", "entropy of formation", "melting point", "net charge", "glycine", "RGD", "glycine fraction", "seq length"]

fs2.add(glycine)
fs3.add(RGD)

#create extra features to calculate frequency of important sequence
res_seq2 = fs(seq)
test_seq2 = fs1(seq)
glycine_seq = fs2(seq)
GFPGER_seq = fs3(seq)
G_count = glycine_seq.columns()
GFPGER_count = GFPGER_seq.columns()
G_count[0]
glycine_fraction = [a/b for a,b in zip(G_count[0],length_seqs)]



#compact multiple single-value features and export the columns
compact_seq = compact(res_seq2)
export_columns = compact_seq.columns()

#create a dataframe
dataframe = pd.DataFrame(export_columns, index = Features, columns = SequenceIds)
dataframe = dataframe.transpose()
dataframe['glycine_fraction'] = glycine_fraction
dataframe['length']= length_seqs
#dataframe.columns(Features2)
dataframe_filtered = dataframe[(dataframe.glycine_fraction > 0.2) & (dataframe.length > 1000)]
leftover_dataframe = dataframe[(dataframe.glycine_fraction < 0.2) & (dataframe.length < 1000)]

#export to excel file
export_excel = dataframe.to_excel (r'C:/Users/innymp/Downloads/collagen5.xlsx', index = True, header=True)
export_excel = dataframe_filtered.to_excel (r'C:/Users/innymp/Downloads/collagen6.xlsx', index = True, header=True)

#performing correlation matrix to reduce the correlated features
corr = dataframe_filtered.corr()
cmap = sns.diverging_palette(h_neg= 250, h_pos = 15,s=75,l=40,n=20, as_cmap=True)
mask = np.triu(np.ones_like(corr, dtype=bool))
corr_df = dataframe_filtered.corr().abs()

#explore to see which features to be dropped
#mask = np.triu(np.ones_like(corr_df, dtype=bool))
#tri_df = corr_df.mask(mask)
#chart = sns.heatmap(tri_df,mask=mask,center = 0, cmap=cmap, linewidths=1, annot=True, fmt=".2f")
#chart.set_xticklabels(labels = Features2, rotation=28, fontsize = 12)
#chart.set_yticklabels(labels = Features2, fontsize=14)
#plt.title('Correlation Matrix of the Protein Features',fontsize= 20)
#plt.show()
#to_drop = [c for c in tri_df.columns if any(tri_df[c] > 0.95)]
#to_drop2 = ['length']
#reduced_df = dataframe_filtered.drop(to_drop, axis=1)

#apply dimensionality reduction and PCA analysis
pipe = Pipeline([
        ('scaler', RobustScaler()),
        ('reducer', PCA(n_components = 7,whiten=True))])
pc = pipe.fit_transform(reduced_df)

scaler2 = RobustScaler()
std_df = scaler2.fit_transform(reduced_df)
reducer2 = PCA(n_components = 7)
pa = reducer2.fit_transform(std_df)


#add PC1 and PC2 components to the data to see correlation
reduced_df['PC_1'] = pc[:,0]
reduced_df['PC_2'] = pc[:,1]

#pipe2 = Pipeline([
        #('scaler2', StandardScaler()),
        #('reducer2', PCA(n_components = 10, whiten = True))])

#choose how many components based on 90% variance
print(reducer2.explained_variance_ratio_.cumsum())


var = pipe.steps[1][1].explained_variance_ratio_
#plt.plot(var)
#plt.xlabel('Principal component index', fontsize = 16)
#plt.ylabel('Explained variance ratio', fontsize = 16)
#plt.title('Determining Number of Principal Components based on Variance', fontsize = 20)
#plt.show()

#graph the relationship betwen PC1 and PC2


#mergings = linkage()

array = np.array([pc[:,0],pc[:,1], pc[:,2],pc[:,3]])
array_t = array.T

#use TSNE to inspect data
model = TSNE(n_components=3,perplexity = 10, learning_rate=100)

transformed = model.fit_transform(array_t)

PC1 = transformed[:,0]
PC2 = transformed[:,1]
PC3 = transformed[:,2]

#fig = plt.figure(figsize=(16,10)).gca(projection='3d')

#fig.scatter(
        #xs = array[0],
        #ys = array[1],
        #zs = array[2],
        #cmap='tab10')

#fig.set_xlabel('PC1')
#fig.set_ylabel('PC2')
#fig.set_zlabel('PC3')
#plt.show()

#using hierachical clustering for grouping similarity
mergings = linkage(array_t, method= 'complete', metric = 'Euclidean')

dendrogram(mergings,labels = SequenceIds2_list, leaf_rotation = 90, leaf_font_size=5, color_threshold=3.4)

plt.title("Cluster Analysis of Collagen Type 1 Alpha 1 from Various Organisms", fontsize= 20)

#labels = fcluster(mergings,10, criterion='distance')sns.scatterplot(data=dataframe_filtered , x = 'PC_1',y='PC_2')
plt.show()
