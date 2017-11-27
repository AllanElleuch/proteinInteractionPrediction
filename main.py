# pip install biopython
from Bio import SeqIO
# from sklearn.preprocessing import Normalizer
import pandas
import numpy as np
from sklearn.preprocessing import MinMaxScaler
from sklearn.ensemble import RandomForestClassifier
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
startTime = time.time()
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

def pairSeq(listSeq):
    pairsSeq = []
    for i in range(0,len(listSeq),2):
        pairsSeq.append(listSeq[i]+listSeq[i+1])
    return listSeq


def getFeatures(data,train=True):
    YList = []
    dataHydropathy=[]
    datatension=[]
    dataCharge=[]
    listSeq=[]
    for record in data:
        seq =record.seq.tostring()
        listSeq.append(splitString(seq,8))
        hydropathySequence=[]
        tensionSequence=[]
        chargeSequence=0
        seqSize=len(seq)
        for residue in seq: #X for unknow amino acid residue and U for selenocysteine
            if(residue!=  'X' and residue!=  'U'):
                chargeSequence += calculateChargeAminoAcid(7,residue)
                hydropathySequence.append(hydropathy[residue])
                tensionSequence.append(tension[residue]/seqSize)
        dataCharge.append(chargeSequence)
        dataHydropathy.append(hydropathySequence)
        datatension.append(tensionSequence)

    # print("min : " + str(len(min(listSeq))))

    window_size = 16
    half_window = (window_size-1)//2
    # for seq in dataCharge: # liste avec valeur hydropathy => lissage des valeurs
    #     features.append(sum(seq)/len(seq))
    for tensionSeq,hydropathySeq,chargeSeq in zip(datatension,dataHydropathy,dataCharge): # liste avec valeur hydropathy => lissage des valeurs
        features=[]
        features.append(sum(tensionSeq)) # 0.568for 500 pairs
        features.append(chargeSeq)  #0.596 for 1000 pairs
        features.append(max(hydropathySeq)/len(hydropathySeq)) #54 % for 1000 pairs alone
        YList.append(features)

    pairsSeq = []
    for i in range(0,len(listSeq),2):
        pairsSeq.append(listSeq[i]+listSeq[i+1])
    pairs = pairSeq(listSeq)

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


        # YList.append(sum(seq)/len (seq))

    pairs = []
    for i in range(0,len(YList),2):
        pairs.append(YList[i]+YList[i+1])
    X = np.array(pairs).reshape(-1,len(pairs[0])) # formation en pairs .reshape(-1, 1)

    # if train:
    #     tfidf_matrix = vectorizer.fit_transform(pairsSeq)
    #
    # else :
    #     tfidf_matrix = vectorizer.transform(pairsSeq)
    #
    #
    # final = sparse.hstack((tfidf_matrix, X))

    return X


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

clf = RandomForestClassifier(max_depth=200, random_state=0,n_jobs=-1,n_estimators=400,bootstrap=False,max_features='auto',oob_score = False, verbose = 0)

CROSSVALIDATION = False
GRIDVALIDATION = False
DATASET_E = True  ## DATASET E for testing or split dataset A and B for training and testing
dataSetSize = -1
# dataSetSizeTraining = dataSetSize//2 if dataSetSize >=0 else -1

dataBrutePositive = read("Supp-A-prunned.txt",25000 ) # size 73260
dataBruteNegative = read("Supp-B-prunned.txt",25000 ) # size 72960 #TOTAL 146 220

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
    dataMerge , y_true_merge = mergeInterlacing ( dataBrutePositive , dataBruteNegative ,y , y2 )
else:
    if(not DATASET_E): # Split positive and negative set in 2 for training and testing
        y=[1 for x in range(len(dataBrutePositive)//2)]
        y2=[0 for x in range(len(dataBruteNegative)//2)]
        X , Y = mergeInterlacing ( dataBrutePositive , dataBruteNegative ,y , y2 )
        X_train, y_train  = split_list(X)
        X_test, y_test = split_list(Y)
        dataMerge = X_train
        dataBruteTest = y_train
        y_true_merge = X_test
        yTest =y_test
        print("size dataset for training : " + str(len(X_test)) )
        print("size dataset for testing : " + str(len(y_test)) )
    else:
        y=[1 for x in range(len(dataBrutePositive)//2)]
        y2=[0 for x in range(len(dataBruteNegative)//2)]
        dataMerge , y_true_merge = mergeInterlacing ( dataBrutePositive , dataBruteNegative ,y , y2 )

        dataBruteTest = read("Supp-E-prunned.txt",-1 )
        yTest=[1 for x in range(len(dataBruteTest)//4)]+[0 for x in range(len(dataBruteTest)//4)]
        print("Size Test set " + str(len(yTest)))

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


featuresTest = getFeatures(dataMerge)



features3= getFeatures(dataBruteTest,train=False) # SUppose que données sont des pairs bout à bout





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

print("info1")
print("features shapes :" + str(featuresTest.shape))


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

        scores2 = cross_validate(clf,X= featuresTest,y=np.array(y_true_merge),cv=5,scoring=scoring) # return_train_score=False,
        print(scores2)

        print("Cross validate in 5 k mean ")
        print("mean fit time " + str( mean( scores2['fit_time']))   )
        print("mean test accuracy  " + str( mean( scores2['test_Accuracy']))   )
        print("mean test_Recall  " + str( mean( scores2['test_Recall']))   )
        print("mean test_f1  " + str( mean( scores2['test_f1']))   )
        print("mean test aur_roc  " + str( mean( scores2['test_AUC']))   )
        # print(mean(prediction2))
# print(results)
# import matplotlib.pyplot as plt

#
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


if not CROSSVALIDATION:

    nbFeatures = features3.shape[0]


    print(len(yTest))

    clf.fit(featuresTest,y_true_merge )


    prediction2 = clf.predict(features3)

    print("data test : " +str(len(dataBruteTest)))
    #
    # print(prediction2)
    print("SCORE : " + str(clf.score(features3,yTest)))

    report = str(classification_report(yTest, prediction2,target_names= ["negative","positive"]))
    print(report) # support = nombre dans le Y test  mais pas la prediction

    aur = str(roc_auc_score(yTest, prediction2))
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
