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




def getFeatures(data):
    YList = []
    dataHydropathy=[]
    for record in data:
        seq =record.seq
        y=[]

        for residue in seq: # pass X for unknow amino acid residue and U for selenocysteine
                            #Transform into dataHydropathy
            if(residue!=  'X' and residue!=  'U'):
                y.append(hydropathy[residue])
        dataHydropathy.append(y)
    for seq in dataHydropathy: # liste avec valeur hydropathy => lissage des valeurs
        YList.append(sum(seq)/len(seq))
    return np.array(YList).reshape(-1, 1)


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



dataBrutePositive = read("C:/Users/escroc/Documents/projectBioInformatique/Supp-A-prunned.txt",1000 )
dataBruteNegative = read("C:/Users/escroc/Documents/projectBioInformatique/Supp-B-prunned.txt",1000 )

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

prediction = clf.score(features3,yTest)
print(prediction)
