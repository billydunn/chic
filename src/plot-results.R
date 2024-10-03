library(pROC)
library(plotROC)
library(ggplot2)
library(grid)
library(gridExtra)
library(ComplexUpset)
library(dplyr)
library(stringr)
library(caret)

# This is a script which reads in Random Forest models (10 iterations) and summarises them in three ways
# 1 ROC curves 
# 2 Box plots of the AUC for 10 repeats of each (gene specific) model
# 3 Upset plots to summarise variable importance for the representative "median" model
# 4 Confusion matrices of the representative "median" model using the default 50% probability cut off and the more stringent 90% cut off 

# All of the above shall be repeated for both any size CH and large clone CH 

# output files into any and large ch folders
# Need to edit the script to read an argument --large which will set a boolean variable large to true
# This will also edit the outdir (large), the names of the input files (large instead of any)
# And also the numbers used in accessing the list (2 instead of 1)
# We could add a function that if the script was called without the large argument, it calls itself with the large arg
# This would mean the script recursively does both when called without the large arg

source("functions.R")

# Read in the arguments to the script: for this script to run, this must include an input data frame, a list of genes, and the desired name of the output directory
args <- commandArgs(trailingOnly = TRUE)

# Stop if no arguments provided 
if (length(args) < 1) {
  stop("No command line arguments provided. \nScript usage --genes <geneslist> --outdir <outdirname> --age-sex-match (optional) --large (optional)")
}

if (any(!args[str_detect(args, "--")] %in% c("--genes", "--data", "--outdir", "--age-sex-match", "--large"))) {
  stop(paste0("Unrecognised arguments provided: ", 
              paste0(args[str_detect(args, "--") & !args %in% c("--genes", "--outdir")]),
              " \nScript usage --data <dataset> --genes <geneslist> --outdir <outdirname> --age-sex-match (optional) --large (optional)"))
}

if (!length(args) %in% c(6,7,8)) {
  stop("Incorrect number of command line arguments provided. \nScript usage --data <dataset> --genes <geneslist> --outdir <outdirname> --age-sex-match (optional)")
}

if (any(!args[str_detect(args, "--")] %in% c("--data", "--genes", "--outdir", "--age-sex-match", "--large"))) {
  stop("Incorrect command line arguments provided. \nScript usage --data <dataset> --genes <geneslist> --outdir <outdirname> --age-sex-match (optional)")
}

# Print the arguments
cat("Running plotting script with command-line arguments:", paste(args, collapse = ", "), "\n")

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

large.gene.list <- readLines(paste0(outdir_path,'/tmp/large_gene_list.txt'))
  
if (any(args %in% "--large")) {
  large.any <- 2
  results.dir <- "large_clone_CH"
  model.dir <- "large"
  iterations <- c(large.gene.list, "all_CH")
} else {
  large.any <- 1
  results.dir <- "any_size_clone_CH"
  model.dir <- "any"
  iterations <- c(gene.list, "all_CH")
}

if (large.any == 2) {
  positive.class <- "large_clone"
  negative.class <- "no_large_clone"
  clone.col <- "largeclone01"
  remove.class <- "clone"
} else {
  positive.class <- "clone"
  negative.class <- "no_mutation"
  clone.col <- "clone"
  remove.class <- "largeclone01"
}

cat(paste0("Processing and summarising results for models of all CH in addition to gene specific models of the following genes:\n", paste(gene.list, collapse = "\n"), "\n"))

models.list <- vector(mode = 'list', length = length(iterations))

# Read in any CH gene specific models 

for (i in 1:length(iterations)){
  # If it is all CH, different directory
  if (iterations[i] == "all_CH") {
    file.path <- paste0(outdir_path, '/models/all_CH/')
    
    gene.specific.models <- vector(mode = 'list', length = 10)
    for (j in 1:10){
      gene.specific.models[[j]] <- readRDS(paste0(file.path, 'all_CH_', model.dir,'_ch_RF_', j,'.rds'))
    }
    models.list[[i]] <- gene.specific.models
  } else {
    file.path <- paste0(outdir_path,'/models/gene_specific/')
    
    gene.specific.models <- vector(mode = 'list', length = 10)
    for (j in 1:10){
      gene.specific.models[[j]] <- readRDS(paste0(file.path, iterations[i], '/', iterations[i], '_', model.dir, '_ch_RF_', j,'.rds'))
    }
    models.list[[i]] <- gene.specific.models
  }
}

names(models.list) <- iterations

# Read in the pre-processed datasets

data.list <- vector(mode = 'list', length = length(iterations))

for (k in 1:length(iterations)){
  if (iterations[k] == "all_CH") {
    gene.specific.data <- vector(mode = 'list', length = 10)
    for (m in 1:10){
      gene.specific.data[[m]] <-readRDS(paste0(outdir_path, '/datasets/all_CH/all_CH_preprocessed_and_downsampled_datasets_matched_case_control_numbers_', m, '.rds'))
    }
    data.list[[k]] <- gene.specific.data
  } else {
    gene.specific.data <- vector(mode = 'list', length = 10)
    for (m in 1:10){
      # Create the vector of probabilities
      gene.specific.data[[m]] <-readRDS(paste0(outdir_path, '/datasets/gene_specific/', iterations[k], '/', iterations[k], '_preprocessed_and_downsampled_datasets_gene_specific_matched_case_control_numbers_', m, '.rds'))
    }
    data.list[[k]] <- gene.specific.data
  }
}

names(data.list) <- iterations

# Structure of the list is data.list[[GENE][[repeat.num]][[any.or.large]]$test -

prob.list <- vector(mode = 'list', length = length(iterations))
roc.list <- vector(mode = 'list', length = length(iterations))

for (n in 1:length(models.list)){
  gene.specific.probs <- vector(mode = 'list', length = 10)
  for (o in 1:10){
    tmp.dat.test <- data.list[[n]][[o]][[large.any]]$test
    any_ch_fit <- models.list[[n]][[o]]
    # Create the vector of probabilities
    gene.specific.probs[[o]] <- predict(any_ch_fit, tmp.dat.test, type="prob") # this is the rfProbs object
  }
  prob.list[[n]] <- gene.specific.probs
}

auc.df <- data.frame(matrix(ncol = length(iterations), nrow = 10))
auc.plot <- data.frame(matrix(ncol = 2, nrow = 10*length(iterations)))
colnames(auc.df) <- iterations
colnames(auc.plot) <- c('gene','auc')

counter <- 1 

cat("Calculating the AUC for the test set for each model...\n")

for (p in 1:length(roc.list)){
  gene.specific.rocs <- vector(mode = 'list', length = 10)
  
  for (q in 1:10){
    tmp.dat.test <- data.list[[p]][[q]][[large.any]]$test
    gene.specific.rocs[[q]] <- roc(tmp.dat.test$clone, prob.list[[p]][[q]][,positive.class], levels = c(negative.class, positive.class)) # Need to change levels to reflect large clone too define as variable
    auc.df[q, p] <- auc(gene.specific.rocs[[q]])
    auc.plot[counter, 1] <- iterations[p]
    auc.plot[counter, 2] <- auc.df[q, p]
    counter <- counter + 1
  }
  
  roc.list[[p]] <- gene.specific.rocs
}

names(prob.list) <- iterations

names(roc.list) <- iterations

# Iteratively calculate the AUC and clculate the min, max and median AUC for each gene
summary.df <- data.frame(matrix(ncol = length(iterations), nrow = 3))
rownames(summary.df) <- c('Median_AUC', 'Min_AUC', 'Max_AUC')
colnames(summary.df) <- iterations
summary.df[1,] <- apply(auc.df, 2, median)
summary.df[2,] <- apply(auc.df, 2, min)
summary.df[3,] <- apply(auc.df, 2, max)

# Calculate the index of the model that has a AUC closest to the median 

index.of.median.auc <- c()

for (r in 1:length(roc.list)){
  tmp.auc <- auc.df[,r]
  ordered.auc <- sort(tmp.auc)
  index.near.median <- which(tmp.auc == ordered.auc[6]) # We have ten variables, so the median is the average of 5 and 6 - we take index 6 for plotting purposes
  index.of.median.auc <- c(index.of.median.auc, index.near.median)
}

# Remember, structure of the list is data.list$GENENAME[[repeat.num]][[any.or.large]][[1]]$test - the [[1]] represents the actual data, the [[2]] and [[3]] are means and SDs used in scaling/centering

# For each element of roc list s, create a plot and make a vector which gives index.of.median.auc[s] as color red and rest as color blue

grob.list <- vector(mode = 'list', length = length(iterations))

get_roc_plot_df <- function(roc_list_element, group_number) {
  roc.plot.df <- data.frame(M = roc_list_element$predictor, 
                            D = ifelse(roc_list_element$response %in% c("clone", "large_clone"), 1, 0),
                            grp = rep(group_number, length(roc_list_element$response)))
  return(roc.plot.df)
}

cat("Plotting ROC curves...\n")

for (s in 1:length(roc.list)){
  roc.plot.df <- do.call(rbind,
                         mapply(get_roc_plot_df, roc.list[[s]], seq_along(roc.list[[s]]), SIMPLIFY = FALSE))
  
  grob.list[[s]] <- ggplot(roc.plot.df, aes(d = D, m = M, color = as.factor(grp))) + 
    geom_roc(labels = FALSE, pointsize = 0) + 
    geom_roc(data = roc.plot.df[roc.plot.df$grp == index.of.median.auc[s],], aes(d = D, m = M, color = "red"), labels = FALSE, pointsize = 0) + 
    style_roc(xlab = "1 - Specificity", ylab = "Sensitivity", guide = TRUE) + 
    scale_color_manual(values = c(rep("lightblue", 10), "red")) + 
    geom_abline(slope = 1) +
    theme_classic() + 
    theme(legend.position = 'none', plot.title = element_text(size = 25)) + 
    ggtitle(names(roc.list)[s]) + 
    annotate('text', x = 0.75, y = 0.5, label = paste0('Median AUC: ', round(summary.df[1, s], 3)), size = 8)
}

# Save the plots 

number.of.plots <- ceiling(length(grob.list) / 10) 

dir.create(paste0(outdir_path, "/results/", results.dir, "/roccurves"), recursive = TRUE)
dir.create(paste0(outdir_path, "/results/", results.dir, "/boxplots"), recursive = TRUE)
dir.create(paste0(outdir_path, "/results/", results.dir, "/upsetplots"), recursive = TRUE)
dir.create(paste0(outdir_path, "/results/", results.dir, "/confusion"), recursive = TRUE)

first.plot <- 1 

for (u in 1:number.of.plots){
  # Plot up to 10 curves in a single plot 
  last.plot <- NA
  if((first.plot + 9) > length(grob.list)) {
    last.plot <- length(grob.list)
  } else {
    last.plot <- first.plot + 9
  }
  pdf(paste0(outdir_path, "/results/", results.dir, "/roccurves/", "ROC_curves_combined", u, ".pdf"), width=14, height=8, onefile = FALSE)
  do.call(grid.arrange, c(grob.list[first.plot:last.plot], ncol = 2))
  dev.off()
  first.plot <- first.plot + 10
}

# Save the ROC curves individually too
for (v in 1:length(grob.list)){
  pdf(paste0(outdir_path, "/results/", results.dir, "/roccurves/", "ROC_curve_", iterations[v], ".pdf"), width=14, height=8, onefile = FALSE)
  print(grob.list[[v]])
  dev.off()
}

cat("Plotting box plots...\n")

# Save boxplot of all gene specific CH
pdf(paste0(outdir_path, "/results/", results.dir, "/boxplots/", "Boxplot_gene_specific_models.pdf"), width=14, height=8, onefile = FALSE)
ggplot(data = auc.plot[auc.plot$gene != 'all_CH',], aes (x = gene, y = auc, fill = gene)) + 
  geom_boxplot(outlier.shape = NA) + geom_hline(yintercept = 0.5, color = 'red') + ylim(0, 1) +
  geom_jitter() +
  xlab("") + ylab("AUC") +
  theme_bw() + theme(legend.title = element_blank(), axis.text.x = element_text(angle = 45, vjust = 0.5))
dev.off()

# Summarise variable importance with Upset plots 

get.top.4.variables <- function(model) {
  v <- varImp(model)
  topVars <- data.frame(v$importance) %>% arrange(desc(Overall)) %>% rownames()
  return(topVars[1:4])
}

#topvars <- lapply(gene.specific.models, get.top.3.variables, medians = index.of.median.auc)

topvars <- vector(mode = 'list', length = length(iterations)) 

for (z in 1:length(topvars)){
  tmp.models.list <- models.list[iterations[z]]
  median.model <- tmp.models.list[[1]][[index.of.median.auc[z]]]
  topvars[[z]] <- get.top.4.variables(median.model)
}

setnames <- unique(as.character(unlist(topvars)))

topvardf <- data.frame(matrix(nrow = length(topvars), ncol = length(setnames)))
colnames(topvardf) <- setnames
rownames(topvardf) <- names(topvars)

for (i in 1:nrow(topvardf)){
  # Take the corresponding list element in topvars
  variables <- unlist(topvars[i])
  
  # For the colnames present, mark 1, else mark 0
  topvardf[i, colnames(topvardf) %in% variables] <- 1
  topvardf[i, !colnames(topvardf) %in% variables] <- 0
}

topvardf$gene_name <- iterations

cat("Plotting upset plot of top four variables to summarise variable importance...\n")

pdf(paste0(outdir_path, "/results/", results.dir, "/upsetplots/", "variable_importance_gene_specific_median_models_top4.pdf"), width=14, height=8, onefile = FALSE)
upset(
  topvardf[!iterations %in% c('all_CH'),],
  intersect=colnames(topvardf)[!colnames(topvardf) %in% 'gene_name'],
  base_annotations=list(
    'Intersection size'=(
      intersection_size(
        bar_number_threshold=1,
        color='grey9',
        fill='grey80'
      )
      + geom_text(
        mapping=aes(label=gene_name),
        position=position_stack(),
        na.rm=TRUE,
        hjust=1.5,
        size=3,
        angle = 90
      )
      + ggtitle('Overlap between top 4 variables in gene-specific models of CH')
    )
  ),
  width_ratio=0.15,
  height_ratio=1/4
)
dev.off()

# Output confusion matrices to file with 0.5 and 0.9 thresholds 
cat("Writing confusion matrices for the test set...\n")

for (z in 1:length(prob.list)){
  tmp.conf.mat <- confusionMatrix(as.factor(ifelse(prob.list[[z]][[index.of.median.auc[z]]][,positive.class] >= 0.5, positive.class, negative.class)), data.list[[z]][[index.of.median.auc[z]]][[large.any]]$test$clone)$table
  write.table(tmp.conf.mat, paste0(outdir_path, "/results/", results.dir, "/confusion/", iterations[z], "_confusion_matrix_50_prob.txt"))
}

for (z in 1:length(prob.list)){
  tmp.conf.mat <- confusionMatrix(as.factor(ifelse(prob.list[[z]][[index.of.median.auc[z]]][,positive.class] >= 0.9, positive.class, negative.class)), data.list[[z]][[index.of.median.auc[z]]][[large.any]]$test$clone)$table
  write.table(tmp.conf.mat, paste0(outdir_path, "/results/", results.dir, "/confusion/", iterations[z], "_confusion_matrix_90_prob.txt"))
}

# Output confusion matrices based on "real world" CH prevalence (instead of 50% prevalence)
# Based on prevalence estimated from UKB (likely underestimates)

cat("Writing confusion matrices assuming real world prevalence...\n")

for (z in 1:length(iterations)) {
  if(str_detect(args[which(args == "--data") + 1], ".rds")) {
    ukbb <- readRDS(args[which(args == "--data") + 1])
  } else {
    ukbb <- read.table(args[which(args == "--data") + 1], header = TRUE)
  }
  
  # Remove all CH aside from the type of interest (ensuring that it == large clone or any clone accordingly, since symbol will be non NA value irrespective of vaf)  
  if (iterations[z] != "all_CH") {
    ukbb <- ukbb[is.na(ukbb$symbol) | (ukbb$symbol == iterations[z] & !is.na(ukbb$symbol) & ukbb[,clone.col] == positive.class),]
  } 
  
  # Inflate the test set by adding lots of non CH controls that haven't been in train set - to reflect the "prevalence" of CH in UKBB 
  inflated.test <- data.list[[z]][[index.of.median.auc[z]]][[large.any]]$test
  
  # What is the prevalence of gene CH in UKBB
  ukbb[,clone.col] <- as.factor(ifelse(ukbb[,clone.col] == positive.class, positive.class, negative.class))
  #ukbb$largeclone01 <- as.factor(ifelse(ukbb$largeclone01 == "large_clone", "large_clone", "no_large_clone"))
  #table(ukbb$largeclone01)
  n.cases.ukbb <- table(ukbb[,clone.col])[1]
  n.controls.ukbb <- table(ukbb[,clone.col])[2]
  
  # How many cases in our test set
  #table(inflated.test$clone)
  n.cases.test <- table(inflated.test$clone)[1]
  n.controls.test <- table(inflated.test$clone)[2]
  
  # So we need to sample n.clones.in.test/n.clones.in.whole.ukbb * n.no.clones.in.ukbb, minus n.no.clones.already.in.test, and join these on to the test set. Further, we need to make sure we do not include those in the training set for the model we are applying 
  n.controls.to.sample <- ((n.cases.test/n.cases.ukbb) * n.controls.ukbb) - n.controls.test
  
  ukbb.to.sample.from <- ukbb[!((rownames(ukbb) %in% rownames(models.list[[z]][[index.of.median.auc[z]]]$trainingData)) | ukbb[,clone.col] == positive.class | (rownames(ukbb) %in% rownames(inflated.test))),]
  
  # Randomly sampled to get extra cases, then combine with the test set
  
  ukbb.sampled <- ukbb.to.sample.from[sample.int(nrow(ukbb.to.sample.from), n.controls.to.sample), 
                                      !colnames(ukbb.to.sample.from) %in% c("symbol", remove.class, "largeclone015", "largeclone02")]
  
  colnames(ukbb.sampled)[colnames(ukbb.sampled) == clone.col] <- "clone"
  
  # Make sure both datasets have the same column names/number/order 
  ukbb.sampled <- ukbb.sampled[,colnames(inflated.test)]
  inflated.test <- rbind(inflated.test, ukbb.sampled)
  
  tmp.model.median <- models.list[[z]][[index.of.median.auc[z]]]
  
  # Run the predictions on the inflated test set
  inflated.pred <- predict(tmp.model.median, inflated.test)
  
  # Demonstrate how this would change using a higher probability threshold - e.g. 75%, 90%
  inflated.probs <- predict(tmp.model.median, inflated.test, type = "prob")
  
  probs.50 <- as.factor(ifelse(inflated.probs[,positive.class] >= 0.5, positive.class, negative.class))
  tmp.conf.mat <- confusionMatrix(probs.50, inflated.test$clone)$table
  write.table(tmp.conf.mat, paste0(outdir_path, "/results/", results.dir, "/confusion/", iterations[z], "_confusion_matrix_50_prob_population_prevalence.txt"))
  
  probs.90 <- as.factor(ifelse(inflated.probs[,positive.class] >= 0.9, positive.class, negative.class))
  tmp.conf.mat <- confusionMatrix(probs.90, inflated.test$clone)$table
  write.table(tmp.conf.mat, paste0(outdir_path, "/results/", results.dir, "/confusion/", iterations[z], "_confusion_matrix_90_prob_population_prevalence.txt"))
}

cat("Finished summarising results!\n")


