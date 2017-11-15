# pip install biopython
from Bio import SeqIO
from sklearn.preprocessing import Normalizer
import pandas
import numpy as np
from sklearn.preprocessing import MinMaxScaler
from sklearn.ensemble import RandomForestClassifier
from sklearn import svm
import random
hydropathy = { 'A': 1.8,'R':-4.5,'N':-3.5,'D':-3.5,'C': 2.5,
        'Q':-3.5,'E':-3.5,'G':-0.4,'H':-3.2,'I': 4.5,
        'L': 3.8,'K':-3.9,'M': 1.9,'F': 2.8,'P':-1.6,
        'S':-0.8,'T':-0.7,'W':-0.9,'Y':-1.3,'V': 4.2 }


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




# for record in SeqIO.parse("C:/Users/escroc/Documents/projectBioInformatique/fasta20171101.seq", "fasta"):
#     data.append(record)
# data = list(SeqIO.parse("C:/Users/escroc/Documents/projectBioInformatique/fasta20171101.seq", "fasta"))


def getFeatures(data):
    YList = []
    dataHydropathy=[]
    for record in data:
        seq =record.seq
        y=[]

        for residue in seq: #X for unknow amino acid residue and U for selenocysteine

            if(residue!=  'X' and residue!=  'U'):
                y.append(hydropathy[residue])

        dataHydropathy.append(y)


    window_size = 15
    half_window = (window_size-1)//2

    for seq in dataHydropathy: # liste avec valeur hydropathy => lissage des valeurs
        num_residues = len(seq)
        y_data = []


        for i in range(half_window, num_residues-half_window):
            average_value = 0.0
            for j in range(-half_window, half_window+1):
                average_value += seq[i+j]
            y_data.append(average_value / window_size)

        # YList.append([max(y_data)/len(y_data) ,sum(y_data)/len(y_data) ])
        YList.append(random.randint(1, 10))

        # YList.append(max(y_data)/len(y_data) *100)
        # YList.append(sum(y_data)/len(y_data) *100)

    pairs = []
    for i in range(0,len(YList),2):
        # pairs.append([YList[i]+YList[i+1]])
        pairs.append(YList[i])
        pairs.append(YList[i+1])
    # print(pairs)
    X = np.array(pairs).reshape(-1, 2) # formation en pairs .reshape(-1, 1)
    # scaler = MinMaxScaler(feature_range=(0, 1)) # Réduction
    # rescaledX = scaler.fit_transform(X)


    # print(rescaledX)

    return X


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



dataBrutePositive = read("C:/Users/escroc/Documents/projectBioInformatique/Supp-A-prunned.txt",300 )
dataBruteNegative = read("C:/Users/escroc/Documents/projectBioInformatique/Supp-B-prunned.txt",300 )

features = getFeatures(dataBrutePositive)# SUppose que données sont des pairs bout à bout
y=[1 for x in range(len(features))]
features2 = getFeatures(dataBruteNegative) # SUppose que données sont des pairs bout à bout
y2=[0 for x in range(len(features2))]
clf = RandomForestClassifier(max_depth=100, random_state=0,n_jobs=-1,n_estimators=50,max_features=None)
# clf = svm.SVC()


clf.fit(features,y )
clf.fit(features2,y2 )
# print(features2)

dataBruteTest = read("C:/Users/escroc/Documents/projectBioInformatique/Supp-E-prunned.txt",100 )
features3= getFeatures(dataBruteTest) # SUppose que données sont des pairs bout à bout


nbFeatures = len(features3)
yTest=[1 for x in range(nbFeatures)]
# yTest=[1 for x in range(nbFeatures//2)]+[0 for x in range(nbFeatures//2)]
# print(clf.get_params())


prediction = clf.score(features,y)
print(prediction)
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
