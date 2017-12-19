# Predicting Protein-Protein Interactions using Machine Learning


### Introduction : 

Prediction of protein-protein interaction by using physico chemical property of proteins.

### Project organisation : 

Python file :

- main.py : Used for loading Dataset, get features, doing crossvalidation, and also predict an external dataset
- aaindex.py : Used to computre indices from AAindex from the AAindex database
- proteinCharge.py : Used to compute the charge dataset
- pruning.py : Used to make PAN dataset compatible with a fasta parseur

Java file : 
- PPI_Preprocessor.java : Script made to export PPI interaction from the DIP dataset and make give it the same shape as PAN dataset
Dataset :
- Supp-A.txt : Positive dataset PPI from PAN 
- Supp-B.txt : Negative dataset PPI from PAN PPI
- Supp-E.txt : First half positive and second half negative dataset from PAN 
- preprocessed.fa: 3000 positive PPI proccessed from the DIP dataset

