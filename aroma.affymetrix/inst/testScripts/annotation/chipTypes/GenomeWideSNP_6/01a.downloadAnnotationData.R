library("aroma.core");
library("aroma.affymetrix");
verbose <- Arguments$getVerbose(-8, timestamp=TRUE);

path <- system.file("testScripts/R", package="aroma.affymetrix");
pathname <- file.path(path, "downloadUtils.R");
source(pathname);

ar <- AromaRepository(verbose=TRUE);

verbose && enter(verbose, "Downloading annotation data");

chipType <- "GenomeWideSNP_6";
tags <- "Full";
verbose && cat(verbose, "Chip type: ", chipType);

pathname <- downloadCDF(ar, chipType);
verbose && cat(verbose, "CDF: ", pathname);

pathname <- downloadCDF(ar, chipType, tags=tags);
verbose && cat(verbose, "CDF: ", pathname);

# Affymetrix NetAffx CSV file
path <- "analysis/downloads/na32/genotyping";
pathnameS <- file.path(path, "GenomeWideSNP_6.na32.annot.csv");
pathname <- downloadAffymetrixNetAffxCsvFile(pathnameS);
verbose && cat(verbose, "NetAffx CSV: ", pathname);

pathnameS <- file.path(path, "GenomeWideSNP_6.cn.na32.annot.csv");
pathname <- downloadAffymetrixNetAffxCsvFile(pathnameS);
verbose && cat(verbose, "NetAffx CN CSV: ", pathname);

verbose && exit(verbose);
