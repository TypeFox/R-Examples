library("aroma.core");
verbose <- Arguments$getVerbose(-8, timestamp=TRUE);
ar <- AromaRepository(verbose=TRUE);

verbose && enter(verbose, "Downloading annotation data");

chipType <- "CytoScanHD_Array";
verbose && cat(verbose, "Chip type: ", chipType);

pathname <- downloadCDF(ar, chipType);
verbose && cat(verbose, "CDF: ", pathname);

pathname <- downloadACS(ar, chipType, tags=".*");
verbose && cat(verbose, "ACS: ", pathname);

pathname <- downloadUFL(ar, chipType, tags=".*");
verbose && cat(verbose, "UFL: ", pathname);

pathname <- downloadUGP(ar, chipType, tags=".*");
verbose && cat(verbose, "UGP: ", pathname);

verbose && exit(verbose);
