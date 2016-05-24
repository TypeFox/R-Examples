library("aroma.core");
verbose <- Arguments$getVerbose(-8, timestamp=TRUE);
ar <- AromaRepository(verbose=TRUE);

verbose && enter(verbose, "Downloading annotation data");


##########################################################################
# Annotation data:
# Hs_PromPR_v02/
#   Hs_PromPR_v02.cdf
#   Hs_PromPR_v02.acs
#   Hs_PromPR_v02.acm
#   Hs_PromPR_v02,unique.acp
##########################################################################
chipType <- "Hs_PromPR_v02";
verbose && cat(verbose, "Chip type: ", chipType);

pathname <- downloadCDF(ar, chipType);
verbose && cat(verbose, "CDF: ", pathname);

#pathname <- downloadACS(ar, chipType, tags=".*");
pathname <- downloadACS(ar, chipType);
verbose && cat(verbose, "ACS: ", pathname);

#pathname <- downloadACM(ar, chipType, tags=".*");
pathname <- downloadACM(ar, chipType);
verbose && cat(verbose, "ACM: ", pathname);

pathname <- downloadACP(ar, chipType, tags="unique");
verbose && cat(verbose, "ACP: ", pathname);

verbose && exit(verbose);
