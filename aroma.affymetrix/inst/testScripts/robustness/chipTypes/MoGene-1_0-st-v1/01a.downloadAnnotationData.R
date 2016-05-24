library("aroma.core");
verbose <- Arguments$getVerbose(-8, timestamp=TRUE);
ar <- AromaRepository(verbose=TRUE);

verbose && enter(verbose, "Downloading annotation data");

chipType <- "MoGene-1_0-st-v1";
verbose && cat(verbose, "Chip type: ", chipType);

pathname <- downloadCDF(ar, chipType, tags="r3");
verbose && cat(verbose, "CDF: ", pathname);

pathname <- downloadACS(ar, chipType, tags=".*");
verbose && cat(verbose, "ACS: ", pathname);

verbose && exit(verbose);
