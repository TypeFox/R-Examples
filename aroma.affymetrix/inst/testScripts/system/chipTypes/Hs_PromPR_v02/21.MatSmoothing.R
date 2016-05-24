###########################################################################
# Replication test
#
# Description:
# This test verifies that aroma.affymetrix can reproduce the estimates
# of the MAT (Model-based Analysis of Tiling arrays) algorithm.
#
# Author: Mark Robinson and Henrik Bengtsson
# Created: 2008-12-09
# Last modified: 2012-09-02
###########################################################################
library("aroma.affymetrix");
library("limma");  # makeContrast()
verbose <- Arguments$getVerbose(-20, timestamp=TRUE);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Data set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dataSet <- "GSE24546";
tags <- "testset";
chipType <- "Hs_PromPR_v02";
sampleNamesMap <- c(
  GSM605951="Prec1_MeDNA_Input1",
  GSM605952="Prec1_MeDNA_IP2",
  GSM605953="Prec1_MeDNA_IP1"
);

cdf <- AffymetrixCdfFile$byChipType(chipType);
print(cdf);

csR <- AffymetrixCelSet$byName(dataSet, tags=tags, cdf=cdf);
print(csR);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Normalize the data using the MAT model
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
mn <- MatNormalization(csR);
csM <- process(mn, verbose=more(verbose, 3));
print(csM);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Convert data set such that it maps to the "unique" CDF
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
csU <- convertToUnique(csM, verbose=verbose);
print(csU);

# Rename
setFullNamesTranslator(csU, function(names, ...) { sampleNamesMap[names] });


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# MatSmoothing: Design #1
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
sampleNames <- getNames(csU);

design1 <- makeContrasts(Prec1_MeDNA_IP1-Prec1_MeDNA_Input1, levels=sampleNames);
colnames(design1) <- "Prec1_IP1_minus_Input";
print(design1);

ms1 <- MatSmoothing(csU, design=design1, probeWindow=600, tag="singleIP");
csMS1 <- process(ms1, units=NULL, verbose=verbose);
print(csMS1);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# MatSmoothing: Design #2
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
sampleNames <- getNames(csU);

design2 <- makeContrasts(Prec1_MeDNA_IP1 + Prec1_MeDNA_IP2-Prec1_MeDNA_Input1, levels=sampleNames);
colnames(design2) <- "Prec1_IPs_minus_Input";

ms2 <- MatSmoothing(csU, design=design2, probeWindow=800, tag="multipleIP")
csMS2 <- process(ms2, units=NULL,verbose=verbose)
print(csMS2);
