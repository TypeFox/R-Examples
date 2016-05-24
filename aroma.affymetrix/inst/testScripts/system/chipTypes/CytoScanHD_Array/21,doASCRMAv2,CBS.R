library("aroma.affymetrix");
verbose <- Arguments$getVerbose(-50, timestamp=TRUE);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setting up data set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dataSetName <- "Affymetrix-CytoScanHD";
chipType <- "CytoScanHD_Array";

csR <- AffymetrixCelSet$byName(dataSetName, chipType=chipType);
print(csR);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Allele-specific CRMAv2
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dsNList <- doASCRMAv2(csR, verbose=verbose);
print(dsNList);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Segmentation
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
seg <- CbsModel(dsNList$total);
print(seg);

ce <- ChromosomeExplorer(seg);
print(ce);

process(ce, arrays=1, chromosomes=1:2, verbose=verbose);
