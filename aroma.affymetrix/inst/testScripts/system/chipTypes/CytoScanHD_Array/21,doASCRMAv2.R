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
