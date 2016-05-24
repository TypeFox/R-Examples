library("aroma.core");
verbose <- Arguments$getVerbose(-8, timestamp=TRUE);
ar <- AromaRepository(verbose=TRUE);

verbose && enter(verbose, "Downloading annotation data");

chipType <- "HuEx-1_0-st-v2";
verbose && cat(verbose, "Chip type: ", chipType);

pathname <- downloadCDF(ar, chipType);
verbose && cat(verbose, "CDF: ", pathname);

pathname <- downloadACS(ar, chipType, tags=".*");
verbose && cat(verbose, "ACS: ", pathname);

pathname <- downloadCDF(ar, chipType, tags="coreR3,A20071112,EP");
verbose && cat(verbose, "CDF: ", pathname);

pathname <- downloadCDF(ar, chipType, tags="fullR3,A20071112,EP");
verbose && cat(verbose, "CDF: ", pathname);

verbose && exit(verbose);
