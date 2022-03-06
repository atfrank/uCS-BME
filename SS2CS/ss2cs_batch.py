import os
import sys
from tqdm import tqdm
import pickle
import numpy as np
import pandas as pd
import argparse
from sklearn import preprocessing
from sklearn import metrics
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_absolute_error
from extractFeatures import extractCT, get_resname_char

def load_ss2cs_model(nucleus, DIR_PATH):
  ''' load save model '''
  filename = DIR_PATH + '/model/RF_' + nucleus + '.sav'
  model = pickle.load(open(filename, 'rb'))
  return(model)

def processfile(inFile, DIR_PATH):
    
    # initialize    
    rna = "user"
    nuclei = ["C1'", "C2'", "C3'", "C4'", "C5'","C2","C5","C6","C8", "H1'", "H2'", "H3'","H4'", "H5'","H5''","H2","H5","H6","H8", "N1", "H1", "N3", "H3"]

    # featurization
    features = extractCT(inFile, rna)
    features.drop('i_resname_char', axis=1, inplace=True)

    # fit one hot encoder
    train_X = pd.read_csv(DIR_PATH+"/data/train_X_NEW.csv",sep=' ',header=0)
    train_X = train_X.drop(['id','length','resid'],axis = 1)
    enc = preprocessing.OneHotEncoder(sparse = False)
    enc.fit(train_X)
    
    # fit model for each nucleus type
    results = pd.DataFrame([])
    for nucleus in nuclei:
    # one hot encoding testing data
        features_resname = features.drop(['id', 'length', 'resid'],axis=1)
        features_info = features['length']
        features_resname_enc = pd.DataFrame(enc.transform(features_resname))
        features_enc = pd.concat([features_info, features_resname_enc],axis = 1)

        # model prediction
        model = load_ss2cs_model(nucleus, DIR_PATH)
        y_pred = model.predict(features_enc)

        # format prediction
        output_resname = features['i_resname'].apply(lambda x: get_resname_char(x))
        output_resid = features['resid']
        output_nucleus = pd.Series([nucleus]*len(features))
        output_cs = pd.Series(y_pred)
        output_error = pd.Series(["."]*len(features))
        result = pd.concat([output_resname, output_resid, output_nucleus, output_cs, output_error],axis=1)
        results = pd.concat([results, result],ignore_index=True)
    
    return results

def concat(inFile, outFile, number_of_structures, DIR_PATH):
    ctfiles = ['%s_%s.ct'%(inFile, f+1) for f in range(number_of_structures)]
    print(ctfiles)
    print("Current Output:", os.path.abspath(outFile))
    print("Processing Files:", len(ctfiles))    
    names = ['model', 'resid', 'resname', 'nucleus', 'simcs', 'id']

    df = pd.DataFrame([])
    for model,ctfile in enumerate(ctfiles):
        df_cur = processfile(ctfile, DIR_PATH)
        df_cur["model"] = model+1
        df_cur["id"] = "."
        df = pd.concat([df, df_cur])
    df = df.drop(columns = [2]).reset_index(drop = True).rename(columns = {0:"nucleus", 1: "simcs", "i_resname": "resname"})
    df = df[names]  
    return df

def main():
    # configure parser
    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--input", help="input prefix for secondary structures stored as CT files", required = True)
    parser.add_argument("-n","--number_of_structures", help="number of secondary structures", required = True, type = int)
    parser.add_argument("-o","--output", help="path to save output CSV", required = True)
    parser.add_argument("-s","--ss2cs_path", help="path to SS2CS", required = True)

    # parse command line
    a = parser.parse_args()  

    # initialize
    inFile = a.input
    outFile = a.output
    DIR_PATH = a.ss2cs_path   
    number_of_structures = a.number_of_structures

	# compute and collect chemical shifts
    results = concat(inFile, outFile, number_of_structures, DIR_PATH)
    results = results[['model', 'resid', 'resname', 'nucleus', 'simcs', 'id']]                   
    results.to_csv(outFile, sep=',', header=True, index=False)    
    print("Done")
# %%

    
if __name__ == "__main__":
    main()