library(caret)
library(doParallel)
library(stringr)
library(dplyr)

# This is a script to run models on the whole UKBB cleaned and decorrelated dataset
# For the publication, the dataset was first cleaned to exclude extreme outliers in blood indices, those with missing data, or prevalent blood cancer diagnoses
# The script will produce two models: one trained to classify the presence/absence of all CH, and one trained to classify large clone (VAF >10%) CH

source("functions.R")

# Read in the arguments to the script: for this script to run, this must include an input data frame, a list of genes, and the desired name of the output directory
args <- commandArgs(trailingOnly = TRUE)

# Stop if no arguments provided 
if (length(args) < 1) {
  stop("No command line arguments provided. \nScript usage --data <inputfile> --genes <geneslist> --outdir <outdirname> --age-sex-match (optional)")
}

if (any(!args[str_detect(args, "--")] %in% c("--data", "--genes", "--outdir", "--features", "--age-sex-match"))) {
  stop(paste0("Unrecognised arguments provided: ", 
             paste0(args[str_detect(args, "--") & !args %in% c("--data", "--genes", "--outdir", "--features")]),
             " \nScript usage --data <inputfile> --genes <geneslist> --outdir <outdirname>"))
}

if (!length(args) %in% c(8,9)) {
  stop("Incorrect number of command line arguments provided. \nScript usage --data <inputfile> --genes <geneslist> --features <featureslist> --outdir <outdirname> --age-sex-match (optional)")
}

if (any(!args[str_detect(args, "--")] %in% c("--data", "--genes", "--outdir", "--features", "--age-sex-match"))) {
  stop("Incorrect command line arguments provided. \nScript usage --data <inputfile> --genes <geneslist> --outdir <outdirname> --age-sex-match (optional)")
}

# Print the arguments
cat("Running gene-specific downsampling script with command-line arguments:", paste(args, collapse = ", "), "\n")

# Read in the input data frame - expects a .rds object but delimited file will work with the script if it has a header
if (str_detect(args[which(args == "--data") + 1], ".rds")) {
  ukbb.master <- readRDS(file = args[which(args == "--data") + 1])
} else {
  ukbb.master <- read.table(file = args[which(args == "--data") + 1], header = TRUE)
}

# Read in the gene list 
gene.list <- readLines(args[which(args == "--genes") + 1])

# Read in the features list 
features.list <- readLines(args[which(args == "--features") + 1])

# Do some pre-processing of the master file: exclude columns that are not in the core variables or specified features
ukbb.master <- ukbb.master[,colnames(ukbb.master) %in% c("symbol", "clone", "largeclone01", features.list)]

# For each of the specified features, check there are no NAs; if there are, exclude rows with NAs
rows.to.exclude <- apply(ukbb.master[,colnames(ukbb.master) %in% features.list], 1, function(x) any(is.na(x)))
ukbb.master <- ukbb.master[!rows.to.exclude,]

# Assess the number of individuals with large clones and add those with >40 to the tmp large gene list
n.large.clones <- ukbb.master %>% filter(largeclone01 == "large_clone" & symbol %in% gene.list) %>% group_by(symbol) %>% summarise(n=n()) %>% filter(n >= 40)
large.clones.to.model <- as.character(n.large.clones$symbol)

cat(paste0("Using a master file with the following features:\n", paste(features.list, collapse = "\n"), "\n"))
cat(paste0(as.integer(table(rows.to.exclude)[2]), " rows had missing data for the specified features and were excluded.\nMaster dataset has data for ", as.integer(table(rows.to.exclude)[1]), " individuals."))

# Read in the outdir
outdir_path <- args[which(args == "--outdir") + 1]

seeds.list <- vector(mode =  'list', length = (length(gene.list) + 1)) 

for (n in 1:length(seeds.list)){
  set.seed(n) # By passing n as the seed, we get 35 lots of 10 different seeds, but we keep the 25 sets of 10 seeds reproducible
  seeds.list[[n]] <- sample.int(1000, 10)
}

iterations <- c(gene.list, "all_CH")

if (any(args %in% "--age-sex-match")) {
  outdir_path <- paste0(outdir_path, "/matched")
  
  for (i in 1:length(iterations)) {
    if (iterations[i] == "all_CH") {
      dir.create(file.path(outdir_path,'datasets/all_CH'), recursive = TRUE)
      ukbb.data <- ukbb.master
      # Remove the parameters we do not want to have in the model (i.e. remove all the clone columns aside from the one of interest e.g. clone, largeclone01)
      features.to.remove.anych <- c('symbol', 'largeclone01', 'largeclone015', 'largeclone02')
      features.to.remove.largech <- c('symbol', 'clone', 'largeclone015', 'largeclone02')
      
      ukbb.data.anych <- ukbb.data[,!colnames(ukbb.data) %in% features.to.remove.anych]
      ukbb.data.largech <- ukbb.data[,!colnames(ukbb.data) %in% features.to.remove.largech]
      
      colnames(ukbb.data.largech)[colnames(ukbb.data.largech) %in% 'largeclone01'] <- 'clone'
      
      ukbb.data.largech$clone <- factor(ifelse(ukbb.data.largech$clone == 'large_clone', 'large_clone', 'no_large_clone'))
      
      datasets <- list(ukbb.data.anych, ukbb.data.largech)
      
      names(datasets) <- c('anych', 'largech')
      
      for (j in 1:10){
        # There is a collection of 10 separate seeds for each iteration of i 
        downsample_anyCH <- downsample_age_sex_match_datasets(datasets[[1]], datasets[[1]]$clone, 'clone', seeds.list[[i]][j]) # We pass the same seeds every time so we produce 10 different data sets per gene, but reproducibly 
        downsample_largeCH <- downsample_age_sex_match_datasets(datasets[[2]], datasets[[2]]$clone, 'large_clone', seeds.list[[i]][j]) # We pass the same seeds every time so we produce 10 different data sets per gene, but reproducibly
        
        downsamples <- list(downsample_anyCH, downsample_largeCH)
        
        saveRDS(downsamples, file = paste0(outdir_path, '/datasets/all_CH/all_CH_downsampled_datasets_matched_case_control_numbers_', j,'.rds'))
      }
    } else {
      dir.create(file.path(outdir_path,'datasets/gene_specific', iterations[i]), recursive = TRUE)
      ukbb.data <- ukbb.master[is.na(ukbb.master$symbol) | ukbb.master$symbol == iterations[i],] # Take only the gene specific form of CH or no CH
      
      # Remove the parameters we do not want to have in the model (i.e. remove all the clone columns aside from the one of interest e.g. clone, largeclone01)
      features.to.remove.anych <- c('symbol', 'largeclone01', 'largeclone015', 'largeclone02')
      features.to.remove.largech <- c('symbol', 'clone', 'largeclone015', 'largeclone02')
      
      ukbb.data.anych <- ukbb.data[,!colnames(ukbb.data) %in% features.to.remove.anych]
      ukbb.data.largech <- ukbb.data[,!colnames(ukbb.data) %in% features.to.remove.largech]
      
      colnames(ukbb.data.largech)[colnames(ukbb.data.largech) %in% 'largeclone01'] <- 'clone'
      
      ukbb.data.largech$clone <- factor(ifelse(ukbb.data.largech$clone == 'large_clone', 'large_clone', 'no_large_clone'))
      
      datasets <- list(ukbb.data.anych, ukbb.data.largech)
      
      names(datasets) <- c('anych', 'largech')
      
      for (j in 1:10){
        # There is a collection of 10 separate seeds for each iteration of i 
        downsample_anyCH <- downsample_age_sex_match_datasets(datasets[[1]], datasets[[1]]$clone, 'clone', seeds.list[[i]][j]) # We pass the same seeds every time so we produce 10 different data sets per gene, but reproducibly 
        downsample_largeCH <- downsample_age_sex_match_datasets(datasets[[2]], datasets[[2]]$clone, 'large_clone', seeds.list[[i]][j]) # We pass the same seeds every time so we produce 10 different data sets per gene, but reproducibly
        
        downsamples <- list(downsample_anyCH, downsample_largeCH)
        
        saveRDS(downsamples, file = paste0(outdir_path, '/datasets/gene_specific/', iterations[i],'/', iterations[i], '_downsampled_datasets_gene_specific_matched_case_control_numbers_', j,'.rds'))
      }
    }
  }
} else {
  outdir_path <- paste0(outdir_path, "/unmatched")
  
  for (i in 1:length(iterations)) {
    if (iterations[i] == "all_CH") {
      dir.create(file.path(outdir_path,'datasets/all_CH'), recursive = TRUE)
      ukbb.data <- ukbb.master
      # Remove the parameters we do not want to have in the model (i.e. remove all the clone columns aside from the one of interest e.g. clone, largeclone01)
      features.to.remove.anych <- c('symbol', 'largeclone01', 'largeclone015', 'largeclone02')
      features.to.remove.largech <- c('symbol', 'clone', 'largeclone015', 'largeclone02')
      
      ukbb.data.anych <- ukbb.data[,!colnames(ukbb.data) %in% features.to.remove.anych]
      ukbb.data.largech <- ukbb.data[,!colnames(ukbb.data) %in% features.to.remove.largech]
      
      colnames(ukbb.data.largech)[colnames(ukbb.data.largech) %in% 'largeclone01'] <- 'clone'
      
      ukbb.data.largech$clone <- factor(ifelse(ukbb.data.largech$clone == 'large_clone', 'large_clone', 'no_large_clone'))
      
      datasets <- list(ukbb.data.anych, ukbb.data.largech)
      
      names(datasets) <- c('anych', 'largech')
      
      for (j in 1:10){
        # There is a collection of 10 separate seeds for each iteration of i 
        downsample_anyCH <- downsample_datasets(datasets[[1]], datasets[[1]]$clone, 'clone', seeds.list[[i]][j]) # We pass the same seeds every time so we produce 10 different data sets per gene, but reproducibly 
        downsample_largeCH <- downsample_datasets(datasets[[2]], datasets[[2]]$clone, 'large_clone', seeds.list[[i]][j]) # We pass the same seeds every time so we produce 10 different data sets per gene, but reproducibly
        
        downsamples <- list(downsample_anyCH, downsample_largeCH)
        
        saveRDS(downsamples, file = paste0(outdir_path, '/datasets/all_CH/all_CH_downsampled_datasets_matched_case_control_numbers_', j,'.rds'))
      }
    } else {
      dir.create(file.path(outdir_path,'datasets/gene_specific', iterations[i]), recursive = TRUE)
      ukbb.data <- ukbb.master[is.na(ukbb.master$symbol) | ukbb.master$symbol == iterations[i],] # Take only the gene specific form of CH or no CH
      
      # Remove the parameters we do not want to have in the model (i.e. remove all the clone columns aside from the one of interest e.g. clone, largeclone01)
      features.to.remove.anych <- c('symbol', 'largeclone01', 'largeclone015', 'largeclone02')
      features.to.remove.largech <- c('symbol', 'clone', 'largeclone015', 'largeclone02')
      
      ukbb.data.anych <- ukbb.data[,!colnames(ukbb.data) %in% features.to.remove.anych]
      ukbb.data.largech <- ukbb.data[,!colnames(ukbb.data) %in% features.to.remove.largech]
      
      colnames(ukbb.data.largech)[colnames(ukbb.data.largech) %in% 'largeclone01'] <- 'clone'
      
      ukbb.data.largech$clone <- factor(ifelse(ukbb.data.largech$clone == 'large_clone', 'large_clone', 'no_large_clone'))
      
      datasets <- list(ukbb.data.anych, ukbb.data.largech)
      
      names(datasets) <- c('anych', 'largech')
      
      for (j in 1:10){
        # There is a collection of 10 separate seeds for each iteration of i 
        downsample_anyCH <- downsample_datasets(datasets[[1]], datasets[[1]]$clone, 'clone', seeds.list[[i]][j]) # We pass the same seeds every time so we produce 10 different data sets per gene, but reproducibly 
        downsample_largeCH <- downsample_datasets(datasets[[2]], datasets[[2]]$clone, 'large_clone', seeds.list[[i]][j]) # We pass the same seeds every time so we produce 10 different data sets per gene, but reproducibly
        
        downsamples <- list(downsample_anyCH, downsample_largeCH)
        
        saveRDS(downsamples, file = paste0(outdir_path, '/datasets/gene_specific/', iterations[i],'/', iterations[i], '_downsampled_datasets_gene_specific_matched_case_control_numbers_', j,'.rds'))
      }
    }
  }
}

for(k in 1:length(iterations)){
  for (m in 1:10){
    if (iterations[k] == "all_CH") {
      downsamples <- readRDS(paste0(outdir_path, '/datasets/all_CH/all_CH_downsampled_datasets_matched_case_control_numbers_', m, '.rds'))
      preprocessed.datasets <- lapply(downsamples, preprocess_datasets)
      saveRDS(preprocessed.datasets, file = paste0(outdir_path, '/datasets/all_CH/all_CH_preprocessed_and_downsampled_datasets_matched_case_control_numbers_', m, '.rds'))
      cat(paste0('completed ', iterations[k], ' sample number ', m, '\n')) 
    } else {
      downsamples <- readRDS(paste0(outdir_path, '/datasets/gene_specific/', iterations[k], '/', iterations[k], '_downsampled_datasets_gene_specific_matched_case_control_numbers_', m, '.rds'))
      preprocessed.datasets <- lapply(downsamples, preprocess_datasets)
      saveRDS(preprocessed.datasets, file = paste0(outdir_path, '/datasets/gene_specific/', iterations[k], '/', iterations[k], '_preprocessed_and_downsampled_datasets_gene_specific_matched_case_control_numbers_', m, '.rds'))
      cat(paste0('completed ', iterations[k], ' sample number ', m, '\n')) 
    }
  }
}

# Save the list of genes to build large clone models of:
dir.create(file.path(outdir_path,'tmp'), recursive = TRUE)
write.table(large.clones.to.model, col.names = FALSE, row.names = FALSE, quote = FALSE, file = paste0(outdir_path, "/tmp/large_gene_list.txt"))

cat("Success. Completed data downsampling and partitioning.\n")
