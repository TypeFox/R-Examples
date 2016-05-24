library("aroma.core");
verbose <- Arguments$getVerbose(-8, timestamp=TRUE);
ar <- AromaRepository(verbose=TRUE);

verbose && enter(verbose, "Downloading annotation data");

chipType <- "GenomeWideSNP_6";
tags <- "Full";
verbose && cat(verbose, "Chip type: ", chipType);

pathname <- downloadACS(ar, chipType, tags=".*");
verbose && cat(verbose, "ACS: ", pathname);

pathname <- downloadCDF(ar, chipType, tags=tags);
verbose && cat(verbose, "CDF: ", pathname);

pathname <- downloadUGP(ar, chipType, tags=c(tags, "na31,hg19", ".*"));
verbose && cat(verbose, "UGP: ", pathname);

pathname <- downloadUFL(ar, chipType, tags=c(tags, "na31,hg19", ".*"));
verbose && cat(verbose, "UFL: ", pathname);

verbose && exit(verbose);
