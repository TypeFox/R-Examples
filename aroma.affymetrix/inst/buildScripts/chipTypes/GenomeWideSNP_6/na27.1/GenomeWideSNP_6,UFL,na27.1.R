if (interactive()) savehistory();
library("aroma.affymetrix");
log <- Verbose(threshold=-10, timestamp=TRUE);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Settings
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
naVersion <- "27.1";
user <- "HB";
datestamp <- "20090519";

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
      main=sprintf(".na%s", naVersion),
#      cn=sprintf(".cn.na%s", naVersion)
      cn=sprintf(".cn.na%s", as.integer(naVersion))
  );

  for (key in names(tagsList)) {
    tags <- tagsList[[key]];
    pathname <- AffymetrixNetAffxCsvFile$findByChipType(chipType, tags=tags);
    if (isFile(pathname)) {
      csvList[[key]] <- AffymetrixNetAffxCsvFile(pathname);
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
  ## GenomeWideSNP_6.na27.1.annot.csv:  int [1:934968] 334945 334944 ...
  ## GenomeWideSNP_6.cn.na27.annot.csv: int [1:945826] 935622 935777 ...
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Update the file footer
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
if (!exists("srcFileTags", mode="list")) {
  srcFileTags <- list();
  srcFiles <- c(list(cdf), csvList);
  for (kk in seq(along=srcFiles)) {
    srcFile <- srcFiles[[kk]];
    tags <- list(
      filename=getFilename(srcFile), 
      filesize=getFileSize(srcFile), 
      checksum=getChecksum(srcFile)
    );
    srcFileTags[[kk]] <- tags;
  }
  print(srcFileTags);
}

footer <- readFooter(ufl);
footer$createdOn <- format(Sys.time(), "%Y%m%d %H:%M:%S", usetz=TRUE);
footer$createdBy = list(
  fullname = "Henrik Bengtsson", 
  email = sprintf("%s@%s", "henrik.bengtsson", "aroma-project.org")
);
names(srcFileTags) <- sprintf("srcFile%d", seq(along=srcFileTags));
footer$srcFiles <- srcFileTags;
writeFooter(ufl, footer);

print(ufl);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# WHAT'S NEW:
#
# o na26 -> na27
#   - No changes
# o na24 -> na26
#   - All changes are for SNP units.
#   - There are 6 NspI and 14 StyI changes in SNP fragment lengths,
#     which some are only minor.
#   - There are 1108 more SNPs that now have missing values.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
ufl <- AromaUflFile$byChipType("GenomeWideSNP_6,Full", tags="na27.1");
ufl2 <- AromaUflFile$byChipType("GenomeWideSNP_6,Full", tags="na26");
x <- summaryOfUnits(ufl);
x2 <- summaryOfUnits(ufl2);

print(x);
##                  snp    cnp affxSnp other   total
##  enzyme1-only 246080 451191       0     0  697271
##  enzyme2-only 160899      0       0     0  160899
##  both         522472 494615       0     0 1017087
##  missing        2495     20    3022   621    6158
##  total        931946 945826    3022   621 1881415
print(x2);
##                  snp    cnp affxSnp other   total
##  enzyme1-only 246080 451191       0     0  697271
##  enzyme2-only 160899      0       0     0  160899
##  both         522472 494615       0     0 1017087
##  missing        2495     20    3022   621    6158
##  total        931946 945826    3022   621 1881415

print(x-x2);
##               snp cnp affxSnp other total
##  enzyme1-only   0   0       0     0     0
##  enzyme2-only   0   0       0     0     0
##  both           0   0       0     0     0
##  missing        0   0       0     0     0
##  total          0   0       0     0     0

# Differences
for (cc in 1:nbrOfColumns(ufl)) {
  units <- whichVector(ufl[,cc] != ufl2[,cc], na.rm=TRUE);
  if (length(units) > 0) {
    df <- cbind(units, ufl[units,cc], ufl2[units,cc], 
                       ufl[units,cc]-ufl2[units,cc]);
    colnames(df)[ncol(df)] <- "delta";
    print(df);
  }
}

# Differences in NAs
for (cc in 1:nbrOfColumns(ufl)) {
  units <- whichVector(is.na(ufl[,cc]) != is.na(ufl2[,cc]));
  df <- cbind(units, ufl[units,cc], ufl2[units,cc]);
  str(df);
}
