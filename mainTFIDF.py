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

hydropathy = { 'A': 1.8,'R':-4.5,'N':-3.5,'D':-3.5,'C': 2.5,
        'Q':-3.5,'E':-3.5,'G':-0.4,'H':-3.2,'I': 4.5,
        'L': 3.8,'K':-3.9,'M': 1.9,'F': 2.8,'P':-1.6,
        'S':-0.8,'T':-0.7,'W':-0.9,'Y':-1.3,'V': 4.2 }


vectorizer = TfidfVectorizer(sublinear_tf=True, max_df=0.8,min_df=1,
                             stop_words='english')
def spliter(rec):
    x=""+rec.seq.tostring()
    res=""
    k=7
    for i in range(0,len(x),k):
        res=res+x[i:i+k]+" "
    return res

def getFeatures(data,train=True):
    YList = []
    dataHydropathy=[]

    seqs = []
    for record in data:
        seqs.append(spliter(record))

    print(seqs)
    print(type(seqs[0]))
    if train:
        tfidf_matrix = vectorizer.fit_transform(seqs)
    else :
        tfidf_matrix = vectorizer.transform(seqs)

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


    return np.array(YList).reshape(-1, 1),tfidf_matrix
# np.array(YList).reshape(-1, 1)


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



dataBrutePositive = read("C:/Users/escroc/Documents/projectBioInformatique/Supp-A-prunned.txt",100 )
dataBruteNegative = read("C:/Users/escroc/Documents/projectBioInformatique/Supp-B-prunned.txt",100 )

features,tfidf_matrix1 = getFeatures(dataBrutePositive)# SUppose que données sont des pairs bout à bout
y=[1 for x in range(len(features))]
features2,tfidf_matrix2 = getFeatures(dataBruteNegative) # SUppose que données sont des pairs bout à bout
y2=[0 for x in range(len(features2))]
clf = RandomForestClassifier(max_depth=100, random_state=0,n_jobs=-1,n_estimators=50,max_features=None)
# clf = svm.SVC()

dataBruteTest = read("C:/Users/escroc/Documents/projectBioInformatique/Supp-E-prunned.txt",100 )


clf.fit(tfidf_matrix1,y )
clf.fit(tfidf_matrix2,y2 )
# print(features2)
features3,tfidf_matrix3= getFeatures(dataBruteTest,train=False) # SUppose que données sont des pairs bout à bout

nbFeatures = len(features3)
yTest=[1 for x in range(nbFeatures)]
# yTest=[1 for x in range(nbFeatures//2)]+[0 for x in range(nbFeatures//2)]
# print(clf.get_params())
# print( tfidf_matrix1.shape)
# print( tfidf_matrix2.shape)
# print( tfidf_matrix3.shape)
# print(len(yTest))
prediction = clf.score(tfidf_matrix3,yTest)
print(prediction)
