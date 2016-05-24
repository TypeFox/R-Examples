library("aroma.core");
verbose <- Arguments$getVerbose(-8, timestamp=TRUE);
ar <- AromaRepository(verbose=TRUE);

verbose && enter(verbose, "Downloading annotation data");

chipType <- "Mapping250K_Nsp";
verbose && cat(verbose, "Chip type: ", chipType);

pathname <- downloadCDF(ar, chipType);
verbose && cat(verbose, "CDF: ", pathname);

pathname <- downloadUGP(ar, chipType, tags=".*");
verbose && cat(verbose, "UGP: ", pathname);

verbose && exit(verbose);
