library("aroma.core");
verbose <- Arguments$getVerbose(-8, timestamp=TRUE);
ar <- AromaRepository(verbose=TRUE);

verbose && enter(verbose, "Downloading annotation data");

chipType <- "GenomeWideSNP_6";
tags <- "Full";
verbose && cat(verbose, "Chip type: ", chipType);

pathname <- downloadCDF(ar, chipType);
verbose && cat(verbose, "CDF: ", pathname);

pathname <- downloadCDF(ar, chipType, tags=tags);
verbose && cat(verbose, "CDF: ", pathname);

pathname <- downloadACS(ar, chipType, tags=".*");
verbose && cat(verbose, "ACS: ", pathname);

pathname <- downloadUFL(ar, chipType, tags=".*");
verbose && cat(verbose, "UFL: ", pathname);

pathname <- downloadUFL(ar, chipType, tags=c(tags, ".*"));
verbose && cat(verbose, "UFL: ", pathname);

pathname <- downloadUGP(ar, chipType, tags=".*");
verbose && cat(verbose, "UGP: ", pathname);

pathname <- downloadUGP(ar, chipType, tags=c(tags, ".*"));
verbose && cat(verbose, "UGP: ", pathname);

pathname <- downloadUGC(ar, chipType, tags=c(tags, ".*"));
verbose && cat(verbose, "UGC: ", pathname);

verbose && exit(verbose);
