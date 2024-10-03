library(caret)
library(doParallel)
library(pROC)
library(dplyr)

# This is a script to run models on the whole UKBB cleaned and decorrelated dataset
# The script will produce two models: one trained to classify the presence/absence of any size CH, and one trained to classify large clone (VAF >10%) CH
# The script iterates through each CH gene and produces gene-specific and all-CH Random Forests models

source("functions.R")

# Read in the arguments to the script: for this script to run, this must include an input data frame, a list of genes, and the desired name of the output directory
args <- commandArgs(trailingOnly = TRUE)

# Stop if no arguments provided 
if (length(args) < 1) {
  stop("No command line arguments provided. \nScript usage --genes <geneslist> --outdir <outdirname> --model-type <type> --age-sex-match (optional)")
}

if (any(!args[str_detect(args, "--")] %in% c("--genes", "--outdir", "--age-sex-match", "--model-type"))) {
  stop(paste0("Unrecognised arguments provided: ", 
              paste0(args[str_detect(args, "--") & !args %in% c("--genes", "--outdir")]),
              " \nScript usage --genes <geneslist> --outdir <outdirname> --age-sex-match (optional) --model-type <type>"))
}

if (!length(args) %in% c(6,7)) {
  stop("Incorrect number of command line arguments provided. \nScript usage --genes <geneslist> --outdir <outdirname> --age-sex-match (optional) --model-type <type>")
}

if (any(!args[str_detect(args, "--")] %in% c("--genes", "--outdir", "--age-sex-match", "--model-type"))) {
  stop("Incorrect command line arguments provided. \nScript usage --genes <geneslist> --outdir <outdirname> --age-sex-match (optional) --model-type <type>")
}

model.type <- args[which(args == "--model-type") + 1]

if (any(!model.type %in% c("RF", "DT", "XGB"))){
  stop(paste0("Unsupported model type specified: ", model.type))
}

# Define the hyperparameters to tune and model variable, depending on model type
# This will be passed to the model building function 

if(model.type == "DT") {
  hyperparams <- expand.grid(.cp = c(0.001, 0.005, 0.01, 0.15, 0.2, 0.3, 0.4, 0.5))
} else if (model.type == "RF") {
  hyperparams <- expand.grid(.mtry = c(2, 4, 6, 8, 10))
} else {
  
  hyperparams <- expand.grid(
    nrounds = 1000,
    eta = 0.3, # Learning rate
    max_depth = c(2, 4, 6, 8, 10), # Tree depth
    gamma = 0, # Regularisation parameter
    colsample_bytree = 1, # Proportion of columns to each tree
    min_child_weight = 1,
    subsample = 1 # Proportion of training set subsampled
  )
} 
                        
caret.pkg <- case_when(model.type == "DT" ~ "rpart",
                       model.type == "RF" ~ "rf",
                       TRUE ~ "xgbTree")

###### NEED TO UPDATE THE FUNCTIONS.R TO PASS THEM THE ABOVE PARAMS AND THE MODEL TYPE
###### ?The rf model needs to call the argument ntree = 1000 
###### Will then need to update the wrapper script to take an argument, t, for type, and call this
###### Can make RF default if it isn't called (if you want) 

# Print the arguments
cat("Running model building script with command-line arguments:", paste(args, collapse = ", "), "\n")

# Read in the gene list 
gene.list <- readLines(args[which(args == "--genes") + 1])

# Read in the outdir
outdir_path <- args[which(args == "--outdir") + 1]

iterations <- c(gene.list, "all_CH")

if (any(args %in% "--age-sex-match")) {
  outdir_path <- paste0(outdir_path, "/matched")
} else {
  outdir_path <- paste0(outdir_path, "/unmatched")
}

# Read in the large gene list
large.gene.list <- readLines(paste0(outdir_path,'/tmp/large_gene_list.txt'))

cat(paste0("Building ", model.type," models for all CH in addition to gene specific models of the following genes:\n", paste(gene.list, collapse = "\n"), "\n"))

n_cores <- detectCores() - 1
cl <- makeCluster(n_cores, type = 'FORK')
registerDoParallel(cl)

start_time <- Sys.time()

for(i in 1:length(iterations)) {
  
  if (iterations[i] == "all_CH") {
    # Create output directories
    dir.create(file.path(outdir_path,'models/all_CH'), recursive = TRUE)
    
    # Build the models 
    for(j in 1:10){
      # Read in the input preprocessed dataset - the subdir depends on whether we are doing matched or unmatched 
      preprocessed.datasets <- readRDS(paste0(outdir_path, '/datasets/', iterations[i], '/', iterations[i], '_preprocessed_and_downsampled_datasets_matched_case_control_numbers_', j, '.rds'))
      
      any_ch_fit <- build.models.in.parallel(dataset = preprocessed.datasets[[1]]$train,
                                             pkg = caret.pkg,
                                             model_type = model.type,
                                             hyperparams = hyperparams)
      print(paste0(iterations[i], ' model for any CH built number ', j))
      
      # Save the model
      saveRDS(any_ch_fit, file = paste0(outdir_path, "/models/", iterations[i], "/all_CH_any_ch_RF_", j, ".rds"))
      
      large_ch_fit <- build.models.in.parallel(dataset = preprocessed.datasets[[2]]$train, 
                                               pkg = caret.pkg,
                                               model_type = model.type,
                                               hyperparams = hyperparams)
      print(paste0(iterations[i], ' specific model for large CH built number', j))
      
      saveRDS(large_ch_fit, file = paste0(outdir_path, "/models/", iterations[i], "/all_CH_large_ch_RF_", j, ".rds"))
    }
    
  } else {
    # Create output directories 
    dir.create(file.path(outdir_path,'models/gene_specific', iterations[i]), recursive = TRUE)
    
    # Enter the loop 
    for(j in 1:10){
      # Read in the input preprocessed dataset - the subdir depends on whether we are doing matched or unmatched 
      preprocessed.datasets <- readRDS(paste0(outdir_path, '/datasets/gene_specific/', iterations[i], '/', iterations[i], '_preprocessed_and_downsampled_datasets_gene_specific_matched_case_control_numbers_', j, '.rds'))
      
      any_ch_fit <- build.models.in.parallel(dataset = preprocessed.datasets[[1]]$train,
                                             pkg = caret.pkg,
                                             model_type = model.type,
                                             hyperparams = hyperparams)
      print(paste0(iterations[i], ' specific model for any CH built number ', j))
      
      # Save the models
      saveRDS(any_ch_fit, file = paste0(outdir_path, "/models/gene_specific/", iterations[i], "/", iterations[i], "_any_ch_RF_", j, ".rds"))
      
      if (iterations[i] %in% large.gene.list) {
        large_ch_fit <- build.models.in.parallel(dataset = preprocessed.datasets[[2]]$train, 
                                                 pkg = caret.pkg,
                                                 model_type = model.type,
                                                 hyperparams = hyperparams)
        print(paste0(iterations[i], ' model for large CH built number', j))
        
        saveRDS(large_ch_fit, file = paste0(outdir_path, "/models/gene_specific/", iterations[i], "/", iterations[i], "_large_ch_RF_", j, ".rds"))
      }
    }  
  }
}

end_time <- Sys.time()

stopCluster(cl)

cat(paste0('Success. Finished model building script on ', end_time))
