# A script to install packages used by the train-ch-models.sh script
cran.packages <- c("stringr",
                   "caret",
                   "doParallel",
                   "dplyr",
                   "pROC",
                   "plotROC",
                   "ggplot2",
                   "grid",
                   "gridExtra",
                   "ComplexUpset")

new.packages <- cran.packages[!(cran.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
