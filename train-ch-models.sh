#!/bin/bash

# Script written by William Dunn 11/9/23

# This is a wrapper script to call scripts to 1) downsample/partition to produce 10 datasets, 2) train ML models and 3) summarise these models and produce plots

# By default, the script produces results based on unmatched input data (that is, does not perform age and sex matching)

# If the -m argument is provided (=match), it will call all the scripts again, but now performing age and sex matching on the data used to train models

# Currently, due to the inefficient loop used to do age and sex matching, the -m option substantially increased the computational time

# The script assumes that CH is denoted clone, or largeclone01 for large clone, and that age and sex are denoted age_at_recruitment and Sex respectively

# Initialize variables with default values
AGESEXMATCH=false

# Describe usage
usage() {
  echo "Usage: $(basename "$0") [options] argument"
  echo
  echo -e "This is a script written by Dr William Dunn (wd289@cam.ac.uk)\nto generate Random Forest binary classifiers to detect the presence of\nclonal haematopoiesis using blood variables"
  echo
  echo "Options:"
  echo -e " -h\tDisplay this help message"
  echo -e " -d\tInput data, either a table or .rds object, see documentation for necessary format and header names (required)"
  echo -e " -g\tList of gene names for gene specific models, must be a plain text file, each on a new line (required)"
  echo -e " -f\tList of features to build models with, must be a plain text file, each feature on a new line (required)"
  echo -e " -o\tDesired name of output directory (required)"
  echo -e " -m\tUse this option if you wish to output models produced after age- and sex-matching between cases and controls (optional)"
  echo
  echo "Example command:"
  echo "  $(basename "$0") -d ukbbdata.rds -g geneslist.txt -f features.txt -o outdir -m"
  echo
  exit 1
}

# Check arguments
while getopts ":mhd:g:f:o:t:" opt; do
  case "$opt" in
    d)
      DATA="$OPTARG"
      ;;
    g)
      GENES="$OPTARG"
      ;;
    f)
      FEATURES="$OPTARG"
      ;;
    o)
      OUTDIR="$OPTARG"
      ;;
    t)
      TYPE="$OPTARG"
      ;;
    m)
      AGESEXMATCH=true
      ;;
    h)
      usage
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      exit 1
      ;;
  esac
done

# Shift the options out, leaving only non-option arguments
shift "$((OPTIND-1))"

# Check if there are any non-option arguments
if [ $# -ne 0 ]; then
  echo "Error: non option arguments found."
  exit 1
fi

# Check the mandatory arguments are used
if [ ! "$DATA" ] || [ ! "$GENES" ] || [ ! "$FEATURES" ] || [ ! "$OUTDIR" ] || [ ! "$TYPE" ]; then
  echo "Missing mandatory arguments: -g, -d, -f, -t and -o must be provided"
  exit 1
fi

if [ "$AGESEXMATCH" = true ]; then
  echo "Running script with age- and sex-matching, this will increase computational time substantially"
else
  echo "Running script in unmatched mode only"
fi

# Check if datasets prepared by looking for .DONE file, if not run the downsampling script

if [ ! -e $OUTDIR/.downsample.DONE ]
then

  Rscript downsample-and-partition.R --data $DATA --genes $GENES --outdir $OUTDIR --features $FEATURES

  touch $OUTDIR/.downsample.DONE

fi

if [ ! -e $OUTDIR/.training.DONE ]
then

  Rscript build-models-10-repeats.R --genes $GENES --outdir $OUTDIR --model-type $TYPE

  touch $OUTDIR/.training.DONE

fi

if [ ! -e $OUTDIR/.plotting.DONE ]
then

  Rscript plot-results.R --data $DATA --genes $GENES --outdir $OUTDIR
  Rscript plot-results.R --large --data $DATA --genes $GENES --outdir $OUTDIR

  touch $OUTDIR/.plotting.DONE

fi

# Repeat the above if age/sex matching requested, with age/sex matching argument
if [ "$AGESEXMATCH" = true ]; then
  if [ ! -e $OUTDIR/.downsample_matched.DONE ]
  then

  Rscript downsample-and-partition.R --data $DATA --genes $GENES --outdir $OUTDIR --features $FEATURES --age-sex-match

  touch $OUTDIR/.downsample_matched.DONE
  fi

  if [ ! -e $OUTDIR/.training_matched.DONE ]
  then

  Rscript build-models-10-repeats.R --genes $GENES --outdir $OUTDIR --age-sex-match --model-type $TYPE

  touch $OUTDIR/.training_matched.DONE
  fi

  if [ ! -e $OUTDIR/.plotting_matched.DONE ]
  then

  Rscript plot-results.R --data $DATA --genes $GENES --outdir $OUTDIR --age-sex-match
  Rscript plot-results.R --large --data $DATA --genes $GENES --outdir $OUTDIR --age-sex-match

  touch $OUTDIR/.plotting_matched.DONE
  fi
fi
