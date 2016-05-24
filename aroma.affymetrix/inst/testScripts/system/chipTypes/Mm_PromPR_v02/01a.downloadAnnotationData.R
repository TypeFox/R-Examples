library("aroma.core");
verbose <- Arguments$getVerbose(-8, timestamp=TRUE);
ar <- AromaRepository(verbose=TRUE);

verbose && enter(verbose, "Downloading annotation data");

chipType <- "Mm_PromPR_v02";
verbose && cat(verbose, "Chip type: ", chipType);

pathname <- downloadCDF(ar, chipType);
verbose && cat(verbose, "CDF: ", pathname);

verbose && exit(verbose);
