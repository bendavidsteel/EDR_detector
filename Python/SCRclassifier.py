import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.svm import SVC 
from sklearn.metrics import classification_report, confusion_matrix

path = "C:/Projects/Matlab/SCR_labelling/"

eda_data = pd.read_table(path + "data.txt")

eda_targets = pd.read_table(path + "targets.txt")

eda_data_train, eda_data_test, eda_targets_train, eda_targets_test = train_test_split(eda_data, eda_targets, test_size = 0.20)

svclassifier = SVC(kernel='sigmoid')  
svclassifier.fit(eda_data_train, eda_targets_train.values.ravel()) 

eda_targets_pred = svclassifier.predict(eda_data_test) 

print confusion_matrix(eda_targets_test, eda_targets_pred)
print classification_report(eda_targets_test, eda_targets_pred)