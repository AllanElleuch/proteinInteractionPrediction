# pip install biopython
from Bio import SeqIO
from sklearn.preprocessing import Normalizer
import pandas
import numpy as np
from sklearn.preprocessing import MinMaxScaler
from sklearn.ensemble import RandomForestClassifier
from sklearn import svm
import random
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.feature_extraction.text import CountVectorizer
import pandas
hydropathy = { 'A': 1.8,'R':-4.5,'N':-3.5,'D':-3.5,'C': 2.5,
        'Q':-3.5,'E':-3.5,'G':-0.4,'H':-3.2,'I': 4.5,
        'L': 3.8,'K':-3.9,'M': 1.9,'F': 2.8,'P':-1.6,
        'S':-0.8,'T':-0.7,'W':-0.9,'Y':-1.3,'V': 4.2 }
# RMN recurent neural network
#keras =>
#https://keras.io/layers/recurrent/

#Faire tfidf concaténer =>  les deux vecteur sur abscisse
#Permuter les valeurs bonne et mauvaise non col 1 et col2 mettre le plus faible en premier et le deuxième plus grand en premier
#Entrainer dans un sens puis dans l'autre
#ngram
#ngram tf idf bigram
mapHydropathy={
'Isoleucine'	:4.5,
'Valine'	:4.2,
'Leucine'	:3.8,
'Phenylalanine':2.8,
'Cysteine'	:2.5,
'Methionine'	:1.9,
'Alanine'	:1.8,
'Glycine'	:-0.4,
'Threonine'	:-0.7,
'Tryptophan'	:-0.9,
'Serine'	:-0.8,
'Tyrosine'	:-1.3,
'Proline'	:-1.6,
'Histidine'	:-3.2,
'Glutamic' 	:-3.5, #acid
'Glutamine'	:-3.5,
'Aspartic' 	:-3.5, #acid
'Asparagine'	:-3.5,
'Lysine'	:-3.9,
'Arginine'	:-4.5
}

charge = {

'A': 58.2, #ALA ,
'G': 79.3, #GLY
'S': 71.7,  #SER
'T': 59.5, #THR
'L': 45.1, #LEU
'I': 46.3, #ILE
'V': 47.5, #VAL
'N': 66, #ASN
'Q': 59.9, #GLN
'R' : 46.2, #ARG
'H': 59.3, #HIS
'W':  52.5, #TRP
'F': 56.2, # PHE
'Y': 59.6,  #TYR
'E':  55.8, #GLU
'D': 61, #ASP
'K':  47.8 , #LYS
'P': 51.4, #PRO
'C':  55.4 , #CYS
'M':  53.4  #MET

}


# for record in SeqIO.parse("C:/Users/escroc/Documents/projectBioInformatique/fasta20171101.seq", "fasta"):
#     data.append(record)
# data = list(SeqIO.parse("C:/Users/escroc/Documents/projectBioInformatique/fasta20171101.seq", "fasta"))

vectorizer = CountVectorizer(analyzer='char_wb', ngram_range=(5, 5)) #ngram_vectorizer
vectorizer = TfidfVectorizer(sublinear_tf=True, max_df=0.8,min_df=1,
                             stop_words='english')
def splitString(x,i):
    res=""
    k=7
    for i in range(0,len(x),k):
        res=res+x[i:i+k]+" "
    return res

def getFeatures(data,train=True):
    YList = []
    dataHydropathy=[]
    dataCharge=[]
    listSeq=[]
    for record in data:
        seq =record.seq.tostring()
        listSeq.append(splitString(seq,8))
        hydropathySequence=[]
        chargeSequence=[]

        for residue in seq: #X for unknow amino acid residue and U for selenocysteine

            if(residue!=  'X' and residue!=  'U'):
                hydropathySequence.append(hydropathy[residue])
                chargeSequence.append(charge[residue])

        dataHydropathy.append(hydropathySequence)
        dataCharge.append(chargeSequence)


    window_size = 16
    half_window = (window_size-1)//2
    features=[]
    # for seq in dataCharge: # liste avec valeur hydropathy => lissage des valeurs
    #     features.append(sum(seq)/len(seq))
    for charged,hydro in zip(dataCharge,dataHydropathy): # liste avec valeur hydropathy => lissage des valeurs
        features.append(sum(charged)/len(charged))
        # features.append(sum(hydro)/len(hydro))
        YList.append(features)

    pairsSeq = []
    for i in range(0,len(listSeq),2):
        pairsSeq.append(listSeq[i]+listSeq[i+1])
    if train:
        tfidf_matrix = vectorizer.fit_transform(pairsSeq)
    else :
        tfidf_matrix = vectorizer.transform(pairsSeq)
    # for seq in dataHydropathy: # liste avec valeur hydropathy => lissage des valeurs
    #     features.append(sum(seq)/len(seq))

        # num_residues = len(seq)

        # hydro_smooth_data = []


        # for i in range(half_window, num_residues-half_window):
        #     average_value = 0.0
        #     for j in range(-half_window, half_window+1):
        #         average_value += seq[i+j]
        #     hydro_smooth_data.append(average_value / window_size)
        # YList.append(sum(y_data)/len(y_data))

        # YList.append([max(y_data)/len(y_data)] )
        # features.append(max(seq)/len(seq))

            # features.append(max(seq)/len(seq))
        # YList.append(random.randint(1, 10))

        # YList.append(max(y_data))


        # YList.append(sum(seq)/len(seq))

    pairs = []
    for i in range(0,len(YList),2):

        pairs.append(YList[i]+YList[i+1])
        # pairs.append(YList[i])
        # pairs.append(YList[i+1])
    # print(pairs)
    # print(pairs)
    X = np.array(pairs).reshape(-1, 2) # formation en pairs .reshape(-1, 1)

    
    # tfidf_matrix
    # scaler = MinMaxScaler(feature_range=(0, 1)) # Réduction
    # rescaledX = scaler.fit_transform(X)
    # print(X.shape)
    # print(X[0])

    # print(rescaledX)
    print(tfidf_matrix)
    print(vectorizer.get_feature_names())
    return tfidf_matrix


# generator = SeqIO.parse("C:/Users/escroc/Documents/projectBioInformatique/fasta20171101.seq", "fasta")
def read(file,number=-1):
    generator = SeqIO.parse(file, "fasta")
    data = []
    if number == -1:
        for record in generator:
            data.append(record)
    else:
        for i in range(number):
            data.append(next(generator))
    return data


dataSetSize = 500
dataBrutePositive = read("C:/Users/escroc/Documents/projectBioInformatique/Supp-A-prunned.txt",dataSetSize)
dataBruteNegative = read("C:/Users/escroc/Documents/projectBioInformatique/Supp-B-prunned.txt",dataSetSize )
dataMerge = dataBrutePositive + dataBruteNegative
print(dataMerge)
# features = getFeatures(dataBrutePositive)# SUppose que données sont des pairs bout à bout
y=[1 for x in range(len(dataBrutePositive)//2)]
# features2 = getFeatures(dataBruteNegative) # SUppose que données sont des pairs bout à bout
y2=[0 for x in range(len(dataBruteNegative)//2)]
clf = RandomForestClassifier(max_depth=100, random_state=0,n_jobs=-1,n_estimators=50,max_features=None)
# clf = svm.SVC()
featuresTest = getFeatures(dataMerge)

# kernel SVM + kernel RDF
from scipy import *
# merge = vstack((features,features2))
# print(merge)
clf.fit(featuresTest,y+y2 )

# clf.fit(features,y )
# clf.fit(features2,y2 )
# print(features2)

dataBruteTest = read("C:/Users/escroc/Documents/projectBioInformatique/Supp-E-prunned.txt",dataSetSize )
features3= getFeatures(dataBruteTest,train=False) # SUppose que données sont des pairs bout à bout


# print(features3)


nbFeatures = features3.shape[0]
# yTest=[1 for x in range(nbFeatures)]
# print(len([1 for x in range(nbFeatures)]))
if dataSetSize==-1:
    yTest=[1 for x in range(nbFeatures//2)]+[0 for x in range(nbFeatures//2)]
else :
    yTest=[1 for x in range(nbFeatures)]

prediction = clf.score(features3,yTest)
print(prediction)

import matplotlib.pyplot as plt
#résultat full data set  0.537194473964 de computation 320.96s] avec 2 features par prot => seq
#0.529224229543 with just mean
#0.515409139214 with max
# plt.hist(features)
# plt.hist(features2)
# plt.show()


# print(len(dataBrutePositive))
# print(len(dataBruteTest))

# res = 0
# for i in prediction:
#     if(i > 0.5):
#         res+=1
# print(res/len(prediction)*100)
# print(prediction)

# print(features2)
# print(features)
# print(features3)


# RandomForestClassifier(bootstrap=True, class_weight=None, criterion='gini',
#             max_depth=2, max_features='auto', max_leaf_nodes=None,
#             min_impurity_decrease=0.0, min_impurity_split=None,
#             min_samples_leaf=1, min_samples_split=2,
#             min_weight_fraction_leaf=0.0, n_estimators=10, n_jobs=1,
#             oob_score=False, random_state=0, verbose=0, warm_start=False)
# print(len(YList[0]))
# import matplotlib.pyplot as plt
# print(plt.style.available)
# plt.plot([x for x in range(len(YList[0]))],YList[0], 'r', marker = 'o')
# plt.style.use('dark_background')
# plt.show()
