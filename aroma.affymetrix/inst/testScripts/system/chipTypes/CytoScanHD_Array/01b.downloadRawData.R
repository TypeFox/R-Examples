path <- system.file("testScripts/R", package="aroma.affymetrix");
pathname <- file.path(path, "downloadUtils.R");
source(pathname);

verbose && enter(verbose, "Downloading raw data");


##########################################################################
# Data set:
# Affymetrix-CytoScanHD/
#  CytoScanHD_Array/
#    []
#
# Overall design:
#
# URL: http://www.affymetrix.com/estore/browse/products.jsp?productId=prod520004#1_3
##########################################################################
dataSet <- "Affymetrix-CytoScanHD";
chipType <- "CytoScanHD_Array";

ds <- downloadAffymetrixDataSet(dataSet, chipType=chipType, pathname="support/downloads/data/cytoscan_demo_data.zip");
print(ds);
## AffymetrixCelSet:
## Name: Affymetrix-HeartBrain
## Tags:
## Path: rawData/Affymetrix-HeartBrain/CytoScanHD_Array
## Platform: Affymetrix
## Chip type: CytoScanHD_Array
## Number of arrays: 7
## Names: 09-1420_B2_Phase4CustomerPanel_CytoScan_PS_20110228, 11-0810_LC_ONC13B_A6_PoP#2_CytoScan-PS_20110511, 11-0816_LC_ONC134B_B10_PoP#2_CytoScan-PS_20110511, ..., C474_A8_CytoScanHD_LabCorp_BetaTest-1_CB_06022011 [7]
## Time period: 2011-03-01 19:07:35 -- 2011-06-07 16:28:18
## Total file size: 460.90MB
## RAM: 0.01MB


verbose && exit(verbose);
