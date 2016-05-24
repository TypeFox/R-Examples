if (interactive()) savehistory();
library("aroma.affymetrix");
log <- Verbose(threshold=-10, timestamp=TRUE);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# WORKAROUNDS
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# The 'GenomeWideSNP_6.cn.na25.annot.csv' [474,602,222 bytes] contains
# a duplicated column header "Fragment Enzyme Type Length Start Stop",
# which needs to be deleted before running this script.  The new file
# length is 474,602,181 bytes. /HB 2008-04-25


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Settings
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
naVersion <- "25";
user <- "HB";
datestamp <- "20080424";

chipType <- "GenomeWideSNP_6";
cdfTags <- "Full";
nbrOfEnzymes <- 2;


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Setup required annotation files
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
if (!exists("cdf")) {
  cdf <- AffymetrixCdfFile$byChipType(chipType, tags=cdfTags);
  rm(csvList);
}
print(cdf);

if (!exists("csvList", mode="list")) {
  csvList <- list();

  tagsList <- c(
#    main=sprintf(   ".na%s", naVersion),
      cn=sprintf(".cn.na%s", naVersion)
  );

  for (key in names(tagsList)) {
    tags <- tagsList[[key]];
    pathname <- AffymetrixNetAffxCsvFile$findByChipType(chipType, tags=tags);
    if (isFile(pathname)) {
      csv <- AffymetrixNetAffxCsvFile(pathname);
      csvList[[key]] <- csv;
      rm(csv);
    }
    rm(tags);
  }
}
print(csvList);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Import UFL from CSV files
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
tags <- sprintf("na%s,%s%s", naVersion, user, datestamp);
ufl <- NULL;
tryCatch({
  ufl <- AromaUflFile$byChipType(getChipType(cdf), tags=tags);
}, error = function(ex) {})
if (is.null(ufl)) {
  ufl <- AromaUflFile$allocateFromCdf(cdf, tags=tags, nbrOfEnzymes=nbrOfEnzymes);
}
print(ufl);

for (kk in seq(along=csvList)) {
  csv <- csvList[[kk]];
  print(csv);
  units <- importFrom(ufl, csv, verbose=log);
  str(units);
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Statistics
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
print(ufl);

x <- summaryOfUnits(ufl, enzymes=c("NspI", "StyI"));
print(x);
##                 snp    cnp affxSnp other   total
## enzyme1-only 246080 451191       0     0  697271
## enzyme2-only 160899      0       0     0  160899
## both         522472 494615       0     0 1017087
## missing        2495     20    3022   621    6158
## total        931946 945826    3022   621 1881415
