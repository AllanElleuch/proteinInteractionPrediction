# pip install biopython
from Bio import SeqIO
# from sklearn.preprocessing import Normalizer
import pandas
import numpy as np
from sklearn.preprocessing import MinMaxScaler
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import RandomForestRegressor

from sklearn import svm
import random
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.feature_extraction.text import CountVectorizer
import pandas
from scipy import sparse
from proteinCharge import *
from sklearn.feature_extraction.text import HashingVectorizer
# from gensim import  models, similarities
import os
from sklearn.model_selection import cross_val_score

from sklearn.feature_extraction.text import CountVectorizer
from sklearn.feature_extraction.text import TfidfTransformer
from sklearn.linear_model import SGDClassifier
from sklearn.model_selection import GridSearchCV
from sklearn.pipeline import Pipeline
from scipy import *

from sklearn.metrics import roc_auc_score
from sklearn import metrics
from sklearn.model_selection import cross_validate
from sklearn.metrics import make_scorer
from sklearn.metrics import accuracy_score
from sklearn.metrics import recall_score
from sklearn.metrics import classification_report
from sklearn.metrics import f1_score
from sklearn.model_selection import train_test_split
import time
import datetime
from ast import literal_eval
startTime = time.time()
hydropathy = { 'A': 1.8,'R':-4.5,'N':-3.5,'D':-3.5,'C': 2.5,
        'Q':-3.5,'E':-3.5,'G':-0.4,'H':-3.2,'I': 4.5,
        'L': 3.8,'K':-3.9,'M': 1.9,'F': 2.8,'P':-1.6,
        'S':-0.8,'T':-0.7,'W':-0.9,'Y':-1.3,'V': 4.2 }
import AAscript

# def path(path):
#     return  open(os.path.join(os.path.dirname(__file__), path))
# RMN recurent neural network
#keras =>
#https://keras.io/layers/recurrent/
#Faire tfidf concaténer =>les deux vecteur sur abscisse
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

tension = {

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

AASCRIPT=True
# for record in SeqIO.parse("C:/Users/escroc/Documents/projectBioInformatique/fasta20171101.seq", "fasta"):
#     data.append(record)
# data = list(SeqIO.parse("C:/Users/escroc/Documents/projectBioInformatique/fasta20171101.seq", "fasta"))

# vectorizer = CountVectorizer(analyzer='char_wb', ngram_range=(5, 5)) #ngram_vectorizer
# preVectorizer = vectorizer = HashingVectorizer(input='content',n_features=500,
#                                        stop_words='english',
#                                        alternate_sign=False, norm='l2',
#                                        binary=False)
tfidf = []

pip = Pipeline([
# ('vect', HashingVectorizer(n_features=3000,ngram_range=(1,5))),
('vect', HashingVectorizer(n_features=3000,ngram_range=(1,2))),
# ('vect', CountVectorizer()),
('tfidf', TfidfTransformer( use_idf=True, smooth_idf=False, sublinear_tf=False)),
# ('clf',TfidfVectorizer(sublinear_tf=True, max_df=0.8,min_df=1,stop_words='english',max_features=500))
# ('tfidf', TfidfTransformer(norm='l2', use_idf=True, smooth_idf=True, sublinear_tf=False)),
])
parameters = {
}
vectorizer=pip
# vectorizer = TfidfVectorizer(sublinear_tf=True, max_df=0.8,min_df=1,stop_words='english',max_features=150) # 964 pour 50 |958  pour 100
# 0.783209351753 avec TfidfVectorizer 0.783209351753% accuracy with 40 000 training sequences and 1800 test sequences
#0.835812964931 avec hashing vectorizer mais prend 2048.301s n_features=3000,ngram_range=(1,5)
#0.768331562168 hashing 524.963s] n_features=500,ngram_range=(1,2)
# 0.774707757705 hashing  749.803s n_features=750,ngram_range=(1,3))
#0.817747077577 hashing 764.563s n_features=1500,ngram_range=(1,2))

def splitString(x,i):
    res=""
    k=7
    for i in range(0,len(x),k):
        res=res+x[i:i+k]+" "
    return res
from operator import add
def pairSeq(listSeq):
    pairsSeq = []
    for i in range(0,len(listSeq),2):
        pairsSeq.append(listSeq[i]+listSeq[i+1])
    return listSeq

def pairSeqMap(listSeq):
    pairsSeq = []
    for i in range(0,len(listSeq),2):
        pairsSeq.append(np.multiply(listSeq[i], listSeq[i+1]) )
    return pairsSeq


def getFeatures(data,train=True):
    YList = []
    dataHydropathy=[]
    datatension=[]
    dataCharge=[]
    dataLength=[]
    hydro_smooth_data = []
    hydro_smooth_pos_data = []
    hydro_smooth_neg_data = []
    listSeq=[]
    protNumber=0
    dataAAINDEX = []
    AAindexList=[]

    for record in data:
        protNumber+=1
        if(protNumber % 50 == 0):
            print("Protein number : " + str(protNumber))
        seq =record.seq.tostring()
        listSeq.append(splitString(seq,8))
        hydropathySequence=[]
        tensionSequence=[]
        chargeSequence=0

        seqSize=len(seq)
        dataLength.append(seqSize)
        # seq2=""
        for residue in seq: #X for unknow amino acid residue and U for selenocysteine
            if(residue!=  'X' and residue!=  'U'):
                # seq2+=residue
                chargeSequence += calculateChargeAminoAcid(7,residue)
                # if(hydropathy[residue] > 0):
                hydropathySequence.append(hydropathy[residue])
                # else:
                #     hydropathySequence.append(0)
                tensionSequence.append(tension[residue])

        if AASCRIPT:
            aaIndex = np.divide(AAscript.getIndexSeq(seq),seqSize)
            data = []
            # print(aaIndex)
            # for item in aaIndex:
            #     data.append(max(item))

            AAindexList.append(aaIndex)


        # np.append(dataCharge,chargeSequence)
        # np.append(dataHydropathy,hydropathySequence)
        # np.append(datatension,tensionSequence)
        # print(datatension)
        # dataCharge.append(calculateCharge(7.35,seq))
        dataHydropathy.append(hydropathySequence)
        datatension.append(tensionSequence)
        dataCharge.append(chargeSequence)

        hydro_pos = []
        hydro_neg = []
        hydro_both = []
        window_size = 21
        half_window = (window_size-1)//2
        num_residues = len(hydropathySequence)
        for i in range(half_window, num_residues-half_window):
            average_value = 0.0
            for j in range(-half_window, half_window+1):
                average_value += hydropathySequence[i+j]
            # if(average_value > 0):
                # hydropathySequence.append(average_value) # 0.0973986680132
            hydro_both.append(average_value)
            if(average_value > 0):
                hydro_pos.append(average_value)
            else:
                hydro_pos.append(0)

            if(average_value < 0):
                hydro_neg.append(average_value)
            else:
                hydro_neg.append(0)
            # else:
            #     hydropathySequence.append(0)

        # hydro_smooth_data.append(sum(hydropathySequence)/num_residues)
        hydro_smooth_data.append(hydro_both)
        hydro_smooth_pos_data.append(hydro_pos)
        hydro_smooth_neg_data.append(hydro_neg)
        # from Bio.SeqUtils.ProtParam import ProteinAnalysis
        # print(seq)
        # analysed_seq = ProteinAnalysis(seq2)
        # hydro_smooth_data.append(analysed_seq.gravy())

    # print("min : " + str(len(min(listSeq))))


    # YList.append(sum(y_data)/len(y_data))
    # for seq in dataCharge: # liste avec valeur hydropathy => lissage des valeurs
    #     features.append(sum(seq)/len(seq)) #REPLACE hydro_smooth_data with dataHydropathy
    for tensionSeq,hydropathySeq,hydropos,hydroneg,chargeSeq, length in zip(datatension,hydro_smooth_data,hydro_smooth_pos_data,hydro_smooth_neg_data,dataCharge,dataLength): # liste avec valeur hydropathy => lissage des valeurs
        features=[] # PARAMETER
        # features.append(sum(tensionSeq)) # 0.568for 500 pairs
        # features.append(chargeSeq/len(tensionSeq))  #0.596 for 1000 pairs
        features.append(sum(tensionSeq)) # 0.568for 500 pairs
        features.append(chargeSeq)  #0.596 for 1000 pairs
        # features.append(sum(hydropathySeq)/len(hydropathySeq)) #54 % for 1000 pairs alone
        # features.append(sum(hydropathySeq)/len(hydropathySeq)) #54 % for 1000 pairs alone
        features.append(sum(hydropathySeq)) #54 % for 1000 pairs alone
        features.append(sum(hydropos)) #54 % for 1000 pairs alone
        features.append(sum(hydroneg)) #54 % for 1000 pairs alone
        features.append(length) #54 % for 1000 pairs alone
        # features.append(sum(hydropathySeq)/len(hydropathySeq)) #54 % for 1000 pairs alone
        YList.append(features)





    pairsSeq = []
    for i in range(0,len(listSeq),2):
        pairsSeq.append(listSeq[i]+listSeq[i+1])

    # for seq in dataHydropathy: # liste avec valeur hydropathy => lissage des valeurs
    #     features.append(sum(seq)/len(seq))

        #





        # YList.append([max(y_data)/len(y_data)] )
        # features.append(max(seq)/len(seq))

            # features.append(max(seq)/len(seq))
        # YList.append(random.randint(1, 10))

        # YList.append(max(y_data))


        # YList.append(sum(seq)/len (seq))

    pairs = []
    # pairsHydroSEQ_A=[]
    # pairsHydroSEQ_B=[]

    for i in range(0,len(YList),2):
        pairs.append(YList[i]+YList[i+1])
        # pairsHydroSEQ_A.append(hydro_smooth_data[i])
        # pairsHydroSEQ_B.append(hydro_smooth_data[i+1])
    print(YList)
    # X = np.array(pairs).reshape(-1,len(pairs[0])) # formation en pairs .reshape(-1, 1)
    X = np.array(pairs)


    if ( TFIDF) :
        if train:
            tfidf_matrix = vectorizer.fit_transform(pairsSeq)

        else :
            tfidf_matrix = vectorizer.transform(pairsSeq)


        final = sparse.hstack((X, tfidf_matrix))
        column = ['tension_A', 'charge_A', 'hydropathy_A','len_A','tension_B', 'charge_B', 'hydropathy_B','len_B','TFIDF', ]
        # X = pandas.DataFrame(final.toarray()) #,  columns=column
        # X = pandas.SparseDataFrame(final.toarray()) #,  columns=column
        # print(type(final))
        # print(final)
        return final
    else:
        column = ['tension_A', 'charge_A', 'hydropathy_both_A','hydropathy_pos_A','hydropathy_neg_A','len_A',
                  'tension_B', 'charge_B', 'hydropathy_both_B','hydropathy_pos_B','hydropathy_neg_B','len_B']
        X2 = pandas.DataFrame(X,  columns=column)
        add_features_dataframe(X2)
        # if(not (PLOTDATA or CROSSVALIDATION) ):
        #     X2 = X2[['hydro_cat','charge_cat','tension_cat']]
        # X = X[['hydro_cat']]
        # X2['hydropathySeq_A'] = pandas.Series(pairsHydroSEQ_A)
        # X2['hydropathySeq_B'] = pandas.Series(pairsHydroSEQ_B)
        # X2['hydropathySeq'] = pandas.Series(dataHydropathy)

        if(AASCRIPT):
            print("SIZE BEFORE " + str(len(AAindexList)))
            AAindexList = pairSeqMap(AAindexList)
            print("SIZE AFTER " + str(len(AAindexList)))

            indices = AAscript.getIndicesUsed()
            dico = {}
            for indice in indices:
                dico[indice]=[]

            for item in AAindexList:
                for i in range(len(indices)):
                    # print(indices[i])
                    # print(dico)
                    # print(dico[indices[i]])
                    # print(AAindexList)
                    dico[indices[i]].append(item[i])
            AAindexList=dico



            # for i in range(len(indices)):
            #     indice = listIndice[i]
            #     valeurIndex = AAindexList[indice]
            featuresList = []
            columns = []
            for key,value in dico.items():
                X2[key]= pandas.Series(value)
            # res=[]
            # for indice in indices:
            #     res.append(dico[indice])
            # AAindexList=res

    # print(AAindexList.shape)



        return X2


# generator = SeqIO.parse("C:/Users/escroc/Documents/projectBioInformatique/fasta20171101.seq", "fasta")
def read(file,number=-1):
    file = os.path.join(os.path.dirname(__file__), file)
    generator = SeqIO.parse(file, "fasta")
    data = []
    if number == -1:
        for record in generator:
            data.append(record)
    else:

        for i in range(number):
            try:
                data.append(next(generator))
            except:
                print("Read error, no line")

    return data

import math
chargeValue = 3

def label_hydro(row, val = chargeValue):
    return row['hydropathy_both_A']*row['hydropathy_both_B']

def label_hydro_pos(row, val = chargeValue):
    # return (row['hydropathy_pos_A']+row['hydropathy_pos_B'])**2/(row['len_B'] +row['len_A'] ) #  TODO
    return (row['hydropathy_pos_A']*row['hydropathy_pos_B'])/(row['len_B'] +row['len_A'] ) #  TODO
    # return (row['hydropathy_pos_A']**2/row['len_A'])  *  (row['hydropathy_pos_B']**2/row['len_B'] ) #  TODO

    # return (row['hydropathy_pos_A']*row['hydropathy_pos_B'])/(row['len_B'] *row['len_A'] ) #  TODO

def label_hydro_neg(row, val = chargeValue):
    # return (row['hydropathy_neg_A']+row['hydropathy_neg_B'])**2/(row['len_B'] +row['len_A'] ) #  TODO
    return (row['hydropathy_neg_A']*row['hydropathy_neg_B'])/(row['len_B'] +row['len_A'] ) #  TODO
    # return (row['hydropathy_neg_A']**2/row['len_A'])  *  (row['hydropathy_neg_B']**2/row['len_B'] ) #  TODO

# aur 0.524866666667 when >0 0.494715555556 when <0 . 0.478262222222 if all



def label_hydro_neg_normalize(row, val = chargeValue):
    # return (row['hydropathy_neg_A']*row['hydropathy_neg_B'])/row['len_B']
    return (row['hydropathy_neg_A']**2/row['len_A'])  *  (row['hydropathy_neg_B']**2/row['len_B'] ) #  TODO
    # return (row['hydropathy_neg_A']row['hydropathy_neg_B'])**2/(row['len_B'] +row['len_A'] ) #  TODO
    # -0.0986538362274 avec carré
    #0.0380956896011 sans carré
    #0.06569190886 sans divisé par len
    # 0.0373495468672 divisé et mettre tout au carré
    # return (row['hydropathy_neg_A']*row['hydropathy_neg_B'])/len(row['length'])
# aur 0.524866666667 when >0 0.494715555556 when <0 . 0.478262222222 if all
def label_hydro_pos_normalize(row, val = chargeValue):
    # return (row['hydropathy_pos_A']/row['len_A'])  *  (row['hydropathy_pos_B']/row['len_B'] ) #  TODO
    # return (row['hydropathy_pos_A']*row['hydropathy_pos_B'])**2/(row['len_B'] +row['len_A'] ) #  TODO
    return (row['hydropathy_pos_A']**2/row['len_A'])  *  (row['hydropathy_pos_B']**2/row['len_B'] ) #  TODO

    # avec caré 0.077923552692
    #avec just div 0.138665348195
    #
    # return (row['hydropathy_pos_A']*row['hydropathy_pos_B'])/len(row['length'])
    #with add 0.5468
    #aur 0.515026666667 with coulon
    #aur 0.509186666667 with div by length of sequence


def label_tension(row, val = chargeValue):
    # return row['tension_A'] * row['tension_B'] #0.0786316120284  # inversé 0.0563729603192
    # return (row['tension_A']+row['tension_B'])**2/(row['len_B'] +row['len_A'] ) #0.126699067382  #dataset inversé 0.0246951352777 avec len
    return (row['tension_A']+row['tension_B'])/(row['len_B'] +row['len_A'] ) #0.126699067382  #dataset inversé 0.0246951352777 avec len 0.042138007757

    #aur 0.529897777778  with *

def label_charge(row, val = chargeValue):
    return (row['charge_A'])*(row['charge_B'])/(row['len_B'] +row['len_A'] ) #0.0814675199656  # dataset inversé 0.0463440044032 sans len et avec 0.069455631278
    # return (row['charge_A']+row['charge_B'])**2/(row['len_B'] +row['len_A'] ) #0.0336219791304
    # return (row['charge_A']+row['charge_B'])**2 #0.0790310621961
    # return (row['charge_A']*row['charge_B'])**2 #0.0490686656137
    # return ( (row['charge_A']/row['len_A']) * (row['charge_B']/row['len_B']) ) # 0.050188599327



def label_hydro_score_negatif(row, val = chargeValue):
    return row['hydropathySeq_A']
    # return literal_eval(row['hydropathySeq_A'])*literal_eval(row['hydropathySeq_B'])
    # return sum(row['hydropathySeq_A'])*sum(row['hydropathySeq_B'])


# def label_hydro_score_positif(row, val = chargeValue):
#     return row['hydropathySeq_A']*row['hydropathySeq_B']

def add_features_dataframe(d):
    import io

    categorieCharge = d.apply (lambda row: label_charge(row,0.1),axis=1)
    categorieHydro = d.apply (lambda row: label_hydro(row,0.1),axis=1)
    categorieHydroPos = d.apply (lambda row: label_hydro_pos(row,0.1),axis=1)
    categorieHydroNeg = d.apply (lambda row: label_hydro_neg(row,0.1),axis=1)
    categorieTension = d.apply (lambda row: label_tension(row,0.1),axis=1)
    # d['hydropathySeq_A']=d['hydropathySeq_A'].apply((lambda x:pandas.eval(x)))#    )eval(d['hydropathySeq_A'])
    # d['hydropathySeq_A']=pandas.eval(d['hydropathySeq_A'],parser='pandas',target=list)

    # apply
    # print(type(d.eval('hydropathySeq_A')[0]))
    # print(type( d['hydropathySeq_A'][0]))

    # print(type(d['hydropathySeq_A'][0]))
    #
    # d['hydropathySeq_A'] =d['hydropathySeq_A'].apply(ast.literal_eval)
    # d['hydropathySeq_B'] =d['hydropathySeq_B'].apply(ast.literal_eval))
    # print(type(d['hydropathySeq_A'][0]))
    # raise
    # categorieHydroNegatif = d.apply (lambda row: label_hydro_score_negatif(row),axis=1)

    # print(categorieHydroNegatif)
    d['charge_cat'] = categorieCharge
    d['hydro_cat'] = categorieHydro
    d['hydro_pos_cat'] = categorieHydroPos
    d['hydro_neg_cat'] = categorieHydroNeg
    d['tension_cat'] = categorieTension
    # d['hydro_neg_cat'] = categorieHydroNegatif
    return d

# clf = RandomForestRegressor(max_depth=None, random_state=0,n_jobs=-1,n_estimators=300,bootstrap=False,max_features='sqrt',oob_score = False, verbose = 0, min_samples_leaf=2,min_samples_split=2)
# clf = RandomForestRegressor(max_depth=None, random_state=0,n_jobs=-1,n_estimators=300,bootstrap=False,max_features='sqrt',oob_score = False, verbose = 0, min_samples_leaf=2,min_samples_split=2)

clf = RandomForestClassifier(max_depth=None, random_state=0,n_jobs=-1,n_estimators=100,bootstrap=True,max_features='auto',oob_score = True, verbose = 0)

CROSSVALIDATION = True
GRIDVALIDATION = False
DATASET_E = True  ## DATASET E for testing or split dataset A and B for training and testing
dataSetSize = -1
TFIDF=False
PLOTDATA = True
useCache = False
SAVE = False
PREDICT=False
# dataSetSizeTraining = dataSetSize//2 if dataSetSize >=0 else -1
featuresList = ['charge_cat','hydro_neg_cat','hydro_pos_cat','tension_cat','KARS160106','MIYS990101','CASG920101','FAUJ880111','FAUJ880112']
# featuresList = ['charge_cat','hydro_neg_cat','hydro_pos_cat','tension_cat']

if not useCache:

    dataBrutePositive = read("Supp-A-prunned.txt",70000 ) # size 73260
    dataBruteNegative = read("Supp-B-prunned.txt",70000 ) # size 72960 #TOTAL 146 220

    print("data dataBrutePositive : " +str(len(dataBrutePositive)))
    print("data dataBruteNegative : " +str(len(dataBruteNegative)))

    def mergeInterlacing(list1,list2, y1,y2):
        X = [None]*(len(list1)+len(list2))
        Y = [None]*(len(y1)+len(y2))
        X[::2] = list1
        X[1::2] = list2
        Y[::2] = y1
        Y[1::2] = y2
        return X,Y
    dataMerge = None
    dataBruteTest= None
    y_true_merge = None
    yTest = None

    def split_list(a_list):
        half = len(a_list)//2
        return a_list[:half], a_list[half:]

    if CROSSVALIDATION: # merge positive and negative for cross validation
        y=[1 for x in range(len(dataBrutePositive)//2)]
        y2=[0 for x in range(len(dataBruteNegative)//2)]
        # dataMerge , y_true_merge = mergeInterlacing ( dataBrutePositive , dataBruteNegative ,y , y2 )
        dataMerge = dataBrutePositive + dataBruteNegative
        y_true_merge = y+y2
    else:
        if(not DATASET_E): # Split positive and negative set in 2 for training and testing
            y=[1 for x in range(len(dataBrutePositive)//2)]
            y2=[0 for x in range(len(dataBruteNegative)//2)]
            # X , Y = mergeInterlacing ( dataBrutePositive , dataBruteNegative ,y , y2 )
            # print("hello")
            # X_train, y_train  = split_list(X)
            # X_test, y_test = split_list(Y)
            # dataMerge = X_train
            # dataBruteTest = y_train
            # y_true_merge = X_test
            #
            # yTest =y_test
            # print("size dataset for training : " + str(len(X_test)) )
            # print("size dataset for testing : " + str(len(y_test)) )
            dataMerge = dataBrutePositive + dataBruteNegative
            y_true_merge = y+y2
        else:
            y=[1 for x in range(len(dataBrutePositive)//2)]
            y2=[0 for x in range(len(dataBruteNegative)//2)]
            # dataMerge , y_true_merge = mergeInterlacing ( dataBrutePositive , dataBruteNegative ,y , y2 )
            # dataMerge , y_true_merge = mergeInterlacing ( dataBrutePositive , dataBruteNegative ,y , y2 )
            dataMerge = dataBrutePositive + dataBruteNegative
            y_true_merge = y+y2
            dataBruteTest = read("Supp-E-prunned.txt",50 ) # TOTO A CHANGER
            yTest=[1 for x in range(len(dataBruteTest)//4)]+[0 for x in range(len(dataBruteTest)//4)]
            print("Size Test set " + str(len(yTest)))

    featuresTrain_set = getFeatures(dataMerge) # training set

    featuresTrain_set['positive_or_negative'] = pandas.Series(y_true_merge)
    print("SAUVEGARDE DU FICHIER train_data_with_features.csv ")
    featuresTrain_set.to_csv('./train_data_with_features.csv', encoding='utf-8',index=False)


    # featuresTest_set= getFeatures(dataBruteTest)
    # featuresTest_set['positive_or_negative'] = pandas.Series(yTest)
    # featuresTest_set.to_csv('./test_data_with_features.csv', encoding='utf-8',index=False)


    # print(len(tab))
    # print(tab[1])
    # for item in tab:
    #     print(tab)


# dataMerge , y_true_merge =  dataBrutePositive + dataBruteNegative ,y +y2
# if dataSetSize==-1:
#
# else :
#     yTest=[1 for x in range(nbFeatures)]

# print("interlacing positive and negative set ")
# print("data merge : " +str(len(dataMerge)))



__location__ = os.path.realpath(
    os.path.join(os.getcwd(), os.path.dirname(__file__)))
def path(filename):
    return  os.path.join(__location__, filename)


import ast

featuresTest_set = pandas.read_csv('test_data_with_features.csv')
add_features_dataframe(featuresTest_set)


# featuresTest_set = pandas.read_csv('train_data_with_features.csv')
# featuresTrain_set = pandas.read_csv('test_data_with_features.csv')


featuresTrain_set = pandas.read_csv('train_data_with_features.csv')
add_features_dataframe(featuresTrain_set)

# if(dataBruteTest):
#     features3= getFeatures(dataBruteTest,train=False) # SUppose que données sont des pairs bout à bout



# raise
    # dataMerge = dataBrutePositive + dataBruteNegative
    # y_true_merge = y+y2
    # dataBruteTest = read("Supp-E-prunned.txt",-1 )
    # yTest
# column = ['tension_A', 'charge_A', 'hydropathy_A','tension_B', 'charge_B', 'hydropathy_B']
# X = pandas.DataFrame(X,  columns=column)
# add_features_dataframe(X)


# raise

# add_features_dataframe(featuresTest)
# featuresTest = featuresTest[['hydro_cat','charge_cat','tension_cat']]

# print(featuresTest)
# featuresTest = featuresTest[['tension_cat']]

# print(featuresTest)





# else
#     dataMerge=dataMerge[:len(dataMerge)/2]
# print("data merge : " +str(len(dataMerge)))
# if(len(dataMerge)%2==0):
    # print("slice : " + str(len(dataMerge)//2))
    # dataMerge=dataMerge[:36480]

# y=[1 for x in range(36480//4)]
# y2=[0 for x in range(36480//4)]



# features = getFeatures(dataBrutePositive)# SUppose que données sont des pairs bout à bout
# features2 = getFeatures(dataBruteNegative) # SUppose que données sont des pairs bout à bout


# print(len(y2))
# param_grid = {
#     # 'n_estimators': [200, 700],
#     # 'max_features': ['auto', 'sqrt', 'log2']
# }
# clf = GridSearchCV(estimator=clf2, param_grid=param_grid, cv= 5)

# clf = svm.SVC()

# print("info1")
# print("features shapes :" + str(featuresTest.shape))


# print("y size : " +  )

# merge = hstack((featuresTest,y+y2))
# np.random.shuffle(merge)
# print("info")
# merge = merge.reshape(featuresTest[0],-1)
# pd = pandas.DataFrame(np.array(featuresTest))
# print(pd)
# print(merge.shape)
# merge = merge.reshape(51,)
# print(merge.shape)


# x=np.hsplit(merge, 1)
# print()
# x2=sum(x[0:-1])

# x=merge[:,[1, 9]]
# print(x)
# print(type(x))
# kernel SVM + kernel RDF
# merge = vstack((features,features2))
# print(merge)


# CROSS VALIDATION TRAINING

# scoring = {'Accuracy': make_scorer(accuracy_score),'Recall': make_scorer(recall_score),'f1':make_scorer(f1_score)}
if CROSSVALIDATION:
    # param_grid = {
    #     'n_estimators': [50,200, 700],
    #     'max_features': ['auto', 'sqrt', 'log2']
    # }
    if GRIDVALIDATION:
        param_grid = {
            'n_estimators': [200,300,400],
            "min_samples_split" : [2,4], #def 2
            "bootstrap": [True, False], #true
              "min_samples_split": [2, 5,10], # 2
               "max_depth": [75, 125,None],
            "max_features": ['auto', 'sqrt', 'log2']
        }
        scoring = {'AUC': 'roc_auc', 'Accuracy': make_scorer(accuracy_score)}


        CV_rfc = GridSearchCV(estimator=clf,
                          param_grid=param_grid,
                          scoring=scoring, cv=5, refit='AUC') #, refit='AUC'
        CV_rfc.fit(featuresTest, y_true_merge)
        results = str(CV_rfc.cv_results_)
        results = "cv_results_ : " + results +"\n"
        bestResults = str(CV_rfc.best_score_ )
        bestResults = "best_score_ : " + bestResults
        bestParam = str(CV_rfc.best_params_)+"\n"
        bestParam = "best_params_ : " + bestParam +"\n"
        print(bestResults)
        print("Best parameters : ")
        print(bestParam)


        save =   open(os.path.join(os.path.dirname(__file__), "crossValidation.txt"), 'w')
        # save.write(results+"\n")
        save.write(bestParam+"\n")
        save.write(bestResults+"\n")
    else:
        # CROSS VALIDATION SCORE
        scoring = {'AUC': 'roc_auc', 'Accuracy': make_scorer(accuracy_score),'Recall': make_scorer(recall_score),'f1': make_scorer(f1_score)}

        # featuresList = ['charge_cat','hydro_neg_cat','hydro_pos_cat','tension_cat']

        scores2 = cross_validate(clf,X= featuresTrain_set[featuresList],y=featuresTrain_set['positive_or_negative'],cv=10,scoring=scoring) # return_train_score=False,
        print(scores2)

        print("Cross validate in 10 k mean ")
        print("mean fit time " + str( mean( scores2['fit_time']))   )
        print("mean test accuracy  " + str( mean( scores2['test_Accuracy']))   )
        print("mean test_Recall  " + str( mean( scores2['test_Recall']))   )
        print("mean test_f1  " + str( mean( scores2['test_f1']))   )
        print("mean test aur_roc  " + str( mean( scores2['test_AUC']))   )


#10 fois

#         mean test accuracy  0.59
# mean test_Recall  0.568
# mean test_f1  0.580417920043
# mean test aur_roc  0.624837777778
# [Finished in 37.89s]
        # print(mean(prediction2))
# print(results)

# ###
# mean fit time 2.49424351454
# mean test accuracy  0.59
# mean test_Recall  0.562666666667
# mean test_f1  0.577898179683
# mean test aur_roc  0.624835555556
# [Finished in 89.886s]

    # if (row['charge_A'])*(row['charge_B']) <=0.2 :
    #     return 1
    #     return (row['charge_A'])*(row['charge_B'])
    # elif  row['charge_A'] <val and row['charge_B'] > val :
    #     return 1
    #     return (row['charge_A'])*(row['charge_B'])
    # else:
    #     return  0

def printCorrelation(indice,coef):

    print("Correlation of "+indice + " : "+ str(coef)   )


if (PLOTDATA):
    import matplotlib.pyplot as plt
    import seaborn as sns
    X_train = featuresTrain_set
    # Y_train = y_true_merge
    # X_train['positive_or_negative'] = pandas.Series(Y_train)
    d=X_train
    # add_features_dataframe(X_train)

    hydroPos_normalize = d.apply(lambda row: label_hydro_pos_normalize(row),axis=1)
    hydroNeg_normalize = d.apply(lambda row: label_hydro_neg_normalize(row),axis=1)

    print("Correlation of hydropathy " + str(np.corrcoef(X_train['positive_or_negative'],X_train['hydro_cat'])[0][1])   )
    print("Correlation of hydropathy negatif " + str(np.corrcoef(X_train['positive_or_negative'],X_train['hydro_neg_cat'])[0][1])   )
    print("Correlation of hydropathy positif " + str(np.corrcoef(X_train['positive_or_negative'],X_train['hydro_pos_cat'])[0][1]))
    print("Correlation of hydropathy normalisé negatif " + str(np.corrcoef(X_train['positive_or_negative'],hydroNeg_normalize)[0][1])   )
    print("Correlation of hydropathy normalisé positif " + str(np.corrcoef(X_train['positive_or_negative'],hydroPos_normalize)[0][1])   )
    print("Correlation between hydropathy normalisé  " + str(np.corrcoef(hydroNeg_normalize,hydroPos_normalize)[0][1])   )
    print("Correlation between hydropathy   " + str(np.corrcoef(X_train['hydro_neg_cat'],X_train['hydro_pos_cat'])[0][1])   )

    print("Correlation of charge " + str(np.corrcoef(X_train['positive_or_negative'],X_train['charge_cat'])[0][1])   )
    print("Correlation of tension " + str(np.corrcoef(X_train['positive_or_negative'],X_train['tension_cat'])[0][1])   )


    # AAindexList = np.array(AAindexList).reshape((len(listIndice),))
    # print(AAindexList.shape)
    listIndice = AAscript.getIndicesUsed()

    for i in range(len(listIndice)):
        indice = listIndice[i]
        valeurIndex = d[indice]
        correlation = np.corrcoef(X_train['positive_or_negative'],valeurIndex)[0][1]
        printCorrelation(indice,correlation)

# Correlation of hydropathy 0.05225263945
# Correlation of hydropathy negatif 0.149693299679
# Correlation of hydropathy positif 0.181894137964
# Correlation of hydropathy normalis� negatif 0.0555663470778
# Correlation of hydropathy normalis� positif 0.121700892231
# Correlation between hydropathy normalis�  -0.00196592230574
# Correlation between hydropathy   -0.0863537413635
# Correlation of charge 0.0708508589547
# Correlation of tension 0.0429862715692
# Correlation of KARS160114 : -0.0307599366609
# Correlation of KARS160106 : -0.0903915995469
# Correlation of MIYS990101 : 0.120855883322
# Correlation of CASG920101 : -0.0469612913214
# Correlation of FAUJ880111 : -0.244421706324
# Correlation of FAUJ880112 : 0.0662432266328

    # print(categorieTension = d.apply (lambda row: label_tension(row,0.1),axis=1)


    # lop = np.corrcoef(X_train['positive_or_negative'],hydroNeg_normalize)[0][1]
    # print(lop)
    # print("Correlation of hydropathy negatif normalize " + str(np.corrcoef(X_train['positive_or_negative'],hydroNeg_normalize) [0][1])   )
    # print("Correlation of hydropathy positif normalize" + str(np.corrcoef(X_train['positive_or_negative'],hydroPos_normalize)  [0][1])   )


    # print("Correlation of Charge-hydrophoby " + str(np.corrcoef(X_train['tension_B'],X_train['charge_B'])[0][1])    )
    # print("Correlation of charge-charge " + str(np.corrcoef(X_train['charge_A'],X_train['charge_B'])[0][1])   )
    # print(np.corrcoef(X_train['positive_or_negative'],X_train['charge_cat']))
    # print(np.corrcoef(X_train['positive_or_negative'],X_train['tension_cat']))


    # pandas.DataFrame.plot.scatter
    # d["chargeRatio"]=(d.charge_A/d.charge_B)
    # categorieCharge = d.apply (lambda row: label_charge(row,0.1),axis=1)
    # categorieHydro = d.apply (lambda row: label_hydro(row,0.1),axis=1)
    # categorieTension = d.apply (lambda row: label_tension(row,0.1),axis=1)
    # #
    # d['charge_cat'] = categorieCharge
    # d['hydro_cat'] = categorieHydro
    # d['tension_cat'] = categorieTension
    # add_features_dataframe(d)
    # print(abs(d["charge_A"].mean()))

    # print(X_train)
    # print(categorieCharge)
    # pandas.DataFrame({'Positive': d.groupby('positive_or_negative').get_group(1).charge_cat,
    #           'Negative':   d.groupby('positive_or_negative').get_group(0).charge_cat}).plot.hist(stacked=True)
    # plt.show()
    # pandas.DataFrame({'Positive': d.groupby('positive_or_negative').get_group(1).charge_hydro,
    #           'Negative':   d.groupby('positive_or_negative').get_group(0).charge_hydro}).plot.hist(stacked=True)
    # plt.show()

    # sns.lmplot(x="charge_B", y="charge_A", hue="positive_or_negative",  data=d);
    # sns.pointplot(x="positive_or_negative", y="hydro_cat",hue="positive_or_negative", data=X_train);
    # plt.show()
    # sns.pointplot(x="positive_or_negative", y="tension_cat",hue="positive_or_negative", data=X_train);
    # plt.show()

    # sns.swarmplot(x="charge_cat", y="charge_hydro", hue="positive_or_negative", data=X_train);
    # pandas.DataFrame({'Positive': d.groupby('positive_or_negative').get_group(1).charge_tension,
    #           'Negative':   d.groupby('positive_or_negative').get_group(0).charge_tension}).plot.hist(stacked=True)
    # plt.show()
    # sns.swarmplot(x="positive_or_negative", y="charge_cat", hue="positive_or_negative", data=X_train);

    # sns.distplot(x='charge_cat',hue='positive_or_negative', data = d);
    # sns.distplot(x='charge_cat', data = d);
    # plt.hist([d["chargeRatio"], d["positive_or_negative"]], color=['r','b'], alpha=0.5)
    # sns.lmplot(x="categorieCharge", y="charge_B", hue="positive_or_negative",  data=d);


    # d.groupby('positive_or_negative').charge_A.hist(stacked=True)
    # pandas.DataFrame({'Positive': d.groupby('positive_or_negative').get_group(1).chargeRatio,
    #           'Negative':   d.groupby('positive_or_negative').get_group(0).chargeRatio}).plot.hist(stacked=True)


    # sns.distplot(x);
    # sns.countplot(x="charge_A", data=X_train, palette="Greens_d");
    # sns.regplot(x="chargeRatio", y="positive_or_negative", data=X_train);
    # sns.lmplot(x="tension_A", y="tension_B", hue="positive_or_negative",  data=d);
    # sns.lmplot(x="charge_A",y="charge_B" ,hue="positive_or_negative",  data=d);
    # sns.regplot(x="charge_A", y="charge_B",data=d, lowess=True)
    # sns.lmplot(x="charge_A", y="charge_B", hue='positive_or_negative', data=d);
    # sns.lmplot(x="charge_cat", y="positive_or_negative", hue='positive_or_negative', data=d);
    # sns.lmplot(x="tension_A", y="tension_B", hue='positive_or_negative', data=d);

    # sns.distplot(d['charge_cat'])
    # sns.pointplot(x="positive_or_negative", y="charge_cat", hue='positive_or_negative', data=d)
    # sns.lmplot(x="charge_cat", y="hydropathy_B", hue='positive_or_negative', data=d);
    # sns.countplot(x="charge_cat", data=X_train,hue='positive_or_negative', palette="Greens_d");
    # sns.regplot(x="positive_or_negative", y="charge_cat", data=X_train);
    print("corre")
    # correlation_matrix = np.corrcoef(d[['charge_cat','positive_or_negative']].values)
    # print(correlation_matrix)
    # print(pandas.DataFrame(d['charge_cat']).corrwith(pandas.DataFrame(d['positive_or_negative'])))
    # sns.barplot(x="positive_or_negative", y="charge_cat", hue="positive_or_negative", data=X_train);
    # make class positive and negative couple


    # df = pandas.DataFrame(X_train)
    # columns = ['a','b','c','a','b','c']

    # plt.hist(featuresTest)
    # plt.hist(features2)
# plt.figure(figsize=(13, 13))
# plt.title("GridSearchCV evaluating using multiple scorers simultaneously",
#           fontsize=16)
#
# plt.xlabel("min_samples_split")
# plt.ylabel("Score")
# plt.grid()
#
# ax = plt.axes()
# ax.set_xlim(0, 402)
# ax.set_ylim(0.73, 1)
#
# # Get the regular numpy array from the MaskedArray
# print(results)
# X_axis = np.array(results['n_estimators'].data, dtype=float)
#
# for scorer, color in zip(sorted(scoring), ['g', 'k']):
#     for sample, style in (('train', '--'), ('test', '-')):
#         sample_score_mean = results['mean_%s_%s' % (sample, scorer)]
#         sample_score_std = results['std_%s_%s' % (sample, scorer)]
#         ax.fill_between(X_axis, sample_score_mean - sample_score_std,
#                         sample_score_mean + sample_score_std,
#                         alpha=0.1 if sample == 'test' else 0, color=color)
#         ax.plot(X_axis, sample_score_mean, style, color=color,
#                 alpha=1 if sample == 'test' else 0.7,
#                 label="%s (%s)" % (scorer, sample))
#
#     best_index = np.nonzero(results['rank_test_%s' % scorer] == 1)[0][0]
#     best_score = results['mean_test_%s' % scorer][best_index]
#
#     # Plot a dotted vertical line at the best score for that scorer marked by x
#     ax.plot([X_axis[best_index], ] * 2, [0, best_score],
#             linestyle='-.', color=color, marker='x', markeredgewidth=3, ms=8)
#
#     # Annotate the best score for that scorer
#     ax.annotate("%0.2f" % best_score,
#                 (X_axis[best_index], best_score + 0.005))
#
# plt.legend(loc="best")
# plt.grid('off')
# plt.show()

# clf.fit(features,y )
# clf.fit(features2,y2 )
# print(features2)




# print("data brute test" + str(len(features3)))




############## DO PREDICTION FOR  TEST SET  + Export results


if not CROSSVALIDATION and not PLOTDATA:

    # nbFeatures = features3.shape[0]
    #
    #
    # print(len(yTest))
    # featuresList = ['charge_cat','hydro_neg_cat','hydro_pos_cat','tension_cat']
    # featuresList = ['hydro_neg_cat','hydro_pos_cat']
    featuresList = ['charge_cat','hydro_neg_cat','hydro_pos_cat','tension_cat']

    clf.fit(featuresTrain_set[featuresList],featuresTrain_set['positive_or_negative'] )


    prediction2 = clf.predict(featuresTest_set[featuresList])

    #
    # print(prediction2)
    # print("SCORE : " + str(clf.score(features3,yTest)))

    report = str(classification_report(featuresTest_set['positive_or_negative'], prediction2,target_names= ["negative","positive"]))
    print(report) # support = nombre dans le Y test  mais pas la prediction

    aur = str(roc_auc_score(featuresTest_set['positive_or_negative'], prediction2))
    aur = "aur score " + aur + "\n"
    print(aur)




    ## EXPORT results
    print( '%s function took %0.3f ms' % ("machine learning",  (time.time() - startTime) ))
    computeTime='%s function took %0.3f ms' % ("machine learning", (time.time() - startTime))
    save =   open(os.path.join(os.path.dirname(__file__), "scorePrediction.txt"), 'w')
    save.write(computeTime)
    save.write(report)
    save.write(aur)









##########################################

# print(np.array(yTest))
# print(prediction2)






# print(prediction)
# print(metrics.matthews_corrcoef(np.array(yTest), prediction2))
# scores2 = cross_val_score(clf, features3, np.array(yTest))
# means = clf.cv_results_['mean_test_score']
# stds = clf.cv_results_['std_test_score']
# print(means)
# print(stds)
# for mean, std, params in zip(means, stds, clf.cv_results_['params']):
#     print("%0.3f (+/-%0.03f) for %r"
#           % (mean, std * 2, params))



# scores2 = cross_validate(clf,X= features3,y=np.array(yTest),cv=2,scoring=scoring ) # return_train_score=False,




# pscore = metrics.accuracy_score(np.array(yTest), prediction)
# score = metrics.f1_score(yTest, prediction)
# print(score)
# print(pscore)
# print(score)
# print( np.mean(scores2))
# print(yTest)
# print("Accuracy: %0.2f (+/- %0.2f)" % (scores2.mean(), scores2.std() * 2))
# print(len(scores2))


#résultat full data set  0.537194473964 de computation 320.96s] avec 2 features par prot => seq
#0.529224229543 with just mean
#0.515409139214 with max



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
