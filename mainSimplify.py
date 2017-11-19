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

# from here http://www.rsc.org/suppdata/lc/b8/b813475h/supplementary_information_prediction_of_surface_tension.pdf
#
charge = {

'A': 58.2, #ALA ,
'G': 79.3, #GLY
'S': 71.7,  #SER
'T': 59.5, #THR
'L': 45.1, #LEU
'I': 46.3, #ILE
'V': 47.5, #VAL
'N': 66, #ASN
'Q': 59.9 #GLN
'R' : 46.2, #ARG
'H': 59.3 #HIS
'W':  52.5 #TRP
'F': 56.2 # PHE
'Y': 59.6  #TYR
'E'  55.8, #GLU
'D' 61, #ASP
'K'  47.8 , #LYS
'P' 51.4, #PRO
'C'  55.4 , #CYS
'M'  53.4  #MET

}

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

feat= features+features2

# features
# features2
# clf.fit(features,y )

clf.fit(list(features)+list(features2),y+y2 )
# print(features2)

dataBruteTest = read("C:/Users/escroc/Documents/projectBioInformatique/Supp-E-prunned.txt",100 )
features3= getFeatures(dataBruteTest) # SUppose que données sont des pairs bout à bout


nbFeatures = len(features3)
yTest=[1 for x in range(nbFeatures)]

prediction = clf.score(features3,yTest)
print(prediction)
