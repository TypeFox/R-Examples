library("aroma.core");
verbose <- Arguments$getVerbose(-8, timestamp=TRUE);
ar <- AromaRepository(verbose=TRUE);

verbose && enter(verbose, "Downloading annotation data");

chipType <- "GenomeWideSNP_5";
tags <- "Full,r2";
verbose && cat(verbose, "Chip type: ", chipType);

pathname <- downloadCDF(ar, chipType, tags=tags);
verbose && cat(verbose, "CDF: ", pathname);

pathname <- downloadACS(ar, chipType, tags=".*");
verbose && cat(verbose, "ACS: ", pathname);

pathname <- downloadUFL(ar, chipType, tags=c(tags, ".*"));
verbose && cat(verbose, "UFL: ", pathname);

pathname <- downloadUGP(ar, chipType, tags=c(tags, ".*"));
verbose && cat(verbose, "UGP: ", pathname);

verbose && exit(verbose);
