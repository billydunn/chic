library(stringr)
library(caret)
library(doParallel)
library(dplyr)
library(pROC)
library(plotROC)
library(ggplot2)
library(grid)
library(gridExtra)
library(ComplexUpset)

# Helper functions for use in scripts 

downsample_datasets <- function(dataset, feature.of.interest, minor.class, n) { # This produces a subsample = in size to the minor class
  # Function takes in the dataset, an integer (n.subsample), a factor vector (feature.of.interest) and character name of the minor class
  set.seed(n)
  
  # Split the dataset to only those without the variable of interest (e.g. any clone, large clone)
  dataset.cases <- dataset[feature.of.interest == minor.class,]
  dataset.controls <- dataset[feature.of.interest != minor.class,]
  
  # Take a random subset of this dataset 
  index.subset <- sample.int(nrow(dataset.controls), nrow(dataset.cases))
  dataset.subset.controls <- dataset.controls[index.subset,]
  
  # Combine the random subset of the data with the cases so we have the same amount of cases (minor class) but fewer controls (major class)
  downsampled.dataset <- rbind(dataset.subset.controls, dataset.cases)
  
  # Rearrange the rows at random 
  downsampled.dataset <- downsampled.dataset[sample.int(nrow(downsampled.dataset), nrow(downsampled.dataset)),]
  return(downsampled.dataset)
}

downsample_age_sex_match_datasets <- function(dataset, feature.of.interest, minor.class, n) { # This produces a subsample = in size to the minor class
  # Function takes in the dataset, an integer (n.subsample), a factor vector (feature.of.interest) and character name of the minor class
  set.seed(n)
  
  # Split the dataset to only those without the variable of interest (e.g. any clone, large clone)
  dataset.cases <- dataset[feature.of.interest == minor.class,]
  dataset.controls <- dataset[feature.of.interest != minor.class,]
  
  downsampled.dataset <- dataset.cases # Initialise df with cases, then we will iteratively add age- and sex-matched controls
  
  # For each case, define a list of potential cases matching their sex and within an age range, then randomly sample from this 
  for (i in 1:nrow(dataset.cases)){
    case.sex <- dataset.cases$sex[i]
    case.age <- dataset.cases$Age_at_recruitment[i]
    potential.controls <- dataset.controls[abs(dataset.controls$Age_at_recruitment - case.age) <= 2 & # within 2 years of age
                                             dataset.controls$sex == case.sex & # same sex
                                             !rownames(dataset.controls) %in% rownames(downsampled.dataset),] # not already in the subsampled dataframe
    # Take a random subset of this dataset 
    index.subset <- sample.int(nrow(potential.controls), 1)
    dataset.subset.controls <- potential.controls[index.subset,]
    # Combine the random subset of the data with the cases so we have the same amount of cases (minor class) but fewer controls (major class)
    downsampled.dataset <- rbind(downsampled.dataset, dataset.subset.controls)
  }
  
  # Rearrange the rows at random 
  downsampled.dataset <- downsampled.dataset[sample.int(nrow(downsampled.dataset), nrow(downsampled.dataset)),]
  return(downsampled.dataset)
}

preprocess_datasets <- function(dataset){
  # Function takes in a decorrelated dataset with NZV removed (so that all datasets definitely have the same set of features)
  # This function used to scale/center but does not do this in the present iteration since we are using tree-based models which do not require this preprocessing
  # At present, the function therefore only partitions the datasets and returns each partition in a list
  set.seed(777)
  
  # Partition data, using only the subset of the data passed to the function, ensuring equal distribution of CH/no CH in test/train 
  trainIndex <- createDataPartition(dataset$clone, 
                                    p = .8, 
                                    list = FALSE, 
                                    times = 1)
  
  datasetTrain <- dataset[trainIndex,]
  datasetTest <- dataset[-trainIndex,]
  
  # Return the partitioned dataset
  trainTest <- list(datasetTrain, datasetTest)
  names(trainTest) <- c('train','test')
  return(trainTest)
}

build.models.in.parallel <- function(dataset, model_type, hyperparams, pkg){
  # Takes a cleaned, pre-processed dataset and calls the train function on it, modelling the binary class "clone"
  set.seed(777)
  
  seeds = vector(mode='list',length=101)
  for (i in 1:100) seeds[[i]] = sample.int(1000,100)
  seeds[[101]] = sample.int(1000,1)
  
  trctrl <- trainControl(method = "repeatedcv", 
                         number = 10, 
                         repeats = 10,
                         seeds = seeds, 
                         savePredictions = "final", # save preds for the optimal tuning parameter
                         classProbs = TRUE,  # class probs in addition to preds
                         summaryFunction = twoClassSummary, # This means tuning will happen based on AUC not accuracy as is default
                         allowParallel = TRUE) # Have turned this on now 
  
  # The below performs grid search to optimise hyperparams with a fixed number of values - if e.g. the features list is much longer than expected, this might not be sufficient 
  if (model_type == "DT") {
    # Train DT
    ch_fit <- train(clone ~., data = dataset, 
                    method = pkg,
                    trControl=trctrl,
                    metric = 'ROC',
                    tuneGrid = hyperparams) 
  } else if (model_type == "RF") {
    # Train RF - needs the ntree param
    ch_fit <- train(clone ~., data = dataset, 
                    method = pkg,
                    trControl=trctrl,
                    metric = 'ROC',
                    ntree = 1000,
                    tuneGrid = hyperparams)
  } else {
    # Train XGB - needs nthread = 1 to allow doParallel 
    ch_fit <- train(clone ~., data = dataset, 
                    method = pkg,
                    trControl=trctrl,
                    metric = 'ROC',
                    tuneGrid = hyperparams,
                    nthread = 1) 
  }
  
  return(ch_fit)
}