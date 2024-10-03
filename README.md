## 🩸 CHIC: Clonal Haematopoiesis Inference from Counts
This repo contains the source code and installation scripts required to run the wrapper script used to generate all tree-based ML models in the paper "Machine learning framework for predicting the presence of high-risk clonal haematopoiesis using complete blood count data: a population-based study of 431,531 UK Biobank participants".

## Get set up
Before trying to run the scripts, run the installation script (contained in the scripts/ directory) to install the relevant packages and dependencies.

## Using the framework 
I recommend using the wrapper script, train-ch-models.sh. Full instructions for using this are available if the script is called with the help (-h) argument: 

```bash
./train-ch-models.sh -h
```

This wrapper script will perform random downsampling and train/test partitioning of input datasets, and will subsequently train binary tree-based classifiers to predict the presence/absence of clonal haematopoiesis. It will also output a summary of the results in the form of ROC curves, boxplots and Upset plots. 

## What input data can I use?
The scripts expect input data in a tabular (e.g. .txt or .csv) or R dataset (data frame saved as a .rds file) format. It is intended for use to train tree-based machine learning classifiers to predict the presence/absence of clonal haematopoiesis, and will do this in a driver-agnostic (treating CH as a single entity) and driver-specific (building driver-gene specific classifiers) manner.

In addition to an input data frame, the script requires a list of gene names in a .txt file, which it will use for building driver gene-specific classifiers. It also requires a list of features in a plain .txt file: to date, the script has focused on using complete blood count (CBC) variables and age/sex as input features, but the wrapper script will accept any feature (numeric or categorical) provided it is specified by the same name in the input data frame. Finally, the type of tree-based classifier to be trained can be specified using the -t argument: one of "DT" (Decision Tree), "RF" (Random Forests) or "XGB" (eXtreme Gradient Boosting). 

In terms of prerequisites, the script currently requires the presence of the following columns in the input data, named exactly as follows: "Age_at_recruitment", "sex", "clone" and "largeclone01" signifying the participant's age, sex, and whether or not they have CH detectable at a VAF of >=2% or >=10% respectively. 

