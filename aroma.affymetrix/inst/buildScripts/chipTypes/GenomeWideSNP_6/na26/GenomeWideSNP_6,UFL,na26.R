if (interactive()) savehistory();
library("aroma.affymetrix");
log <- Verbose(threshold=-10, timestamp=TRUE);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Settings
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
naVersion <- "26";
user <- "HB";
datestamp <- "20080722";

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
      main=sprintf(   ".na%s", naVersion),
      cn=sprintf(".cn.na%s", naVersion)
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
  ## GenomeWideSNP_6.na26.annot.csv:    int [1:934968] 334945 334944 ...
  ## GenomeWideSNP_6.cn.na26.annot.csv: int [1:945826] 935622 935777 ...
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


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Statistics
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
print(ufl);

x <- summaryOfUnits(ufl);
print(x);
## GenomeWideSNP_6,Full,na26,HB20080722:
##                 snp    cnp affxSnp other   total
## enzyme1-only 246080 451191       0     0  697271
## enzyme2-only 160899      0       0     0  160899
## both         522472 494615       0     0 1017087
## missing        2495     20    3022   621    6158
## total        931946 945826    3022   621 1881415


## GenomeWideSNP_6,na26,HB20080821:
##                 snp    cnp affxSnp other   total
## enzyme1-only 240001 451191       0     0  691192
## enzyme2-only 154884      0       0     0  154884
## both         510330 494615       0     0 1004945
## missing        1385     20    3022   621    5048
## total        906600 945826    3022   621 1856069


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# WHAT'S NEW:
#
# o na24 -> na26
#   - All changes are for SNP units.
#   - There are 6 NspI and 14 StyI changes in SNP fragment lengths,
#     which some are only minor.
#   - There are 1108 more SNPs that now have missing values.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
ufl <- AromaUflFile$byChipType("GenomeWideSNP_6,Full", tags="na26");
ufl2 <- AromaUflFile$byChipType("GenomeWideSNP_6,Full", tags="na24");
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
##  enzyme1-only 246340 451191       0     0  697531
##  enzyme2-only 161064      0       0     0  161064
##  both         523155 494615       0     0 1017770
##  missing        1387     20    3022   621    5050
##  total        931946 945826    3022   621 1881415

print(x-x2);
##               snp cnp affxSnp other total
## enzyme1-only -260   0       0     0  -260
## enzyme2-only -165   0       0     0  -165
## both         -683   0       0     0  -683
## missing      1108   0       0     0  1108
## total           0   0       0     0     0

# Differences
for (cc in 1:nbrOfColumns(ufl)) {
  units <- whichVector(ufl[,cc] != ufl2[,cc], na.rm=FALSE);
  df <- cbind(units, ufl[units,cc], ufl2[units,cc], 
                     ufl[units,cc]-ufl2[units,cc]);
  colnames(df)[ncol(df)] <- "delta";
  print(df);
}

##      units length length delta
##  1  159802    611    641   -30
##  2  191853    481   1025  -544
##  3  289265    974    976    -2
##  4  510015   1128   1132    -4
##  5  671142    964    226   738
##  6  879970   1187   1474  -287
##  
##      units length.02 length.02 delta
##  1  159802       385      1015  -630
##  2  191853      1007       953    54
##  3  200432       737       744    -7
##  4  205447       747       742     5
##  5  240336       739       684    55
##  6  351955       534       495    39
##  7  415362       744      1228  -484
##  8  446535       741       744    -3
##  9  671142      1217      1516  -299
##  10 783865      1070       857   213
##  11 786653       279       415  -136
##  12 817336       844       951  -107
##  13 821426       485       637  -152
##  14 879970       740       748    -8

# Differences in NAs
for (cc in 1:nbrOfColumns(ufl)) {
  units <- whichVector(is.na(ufl[,cc]) != is.na(ufl2[,cc]));
  df <- cbind(units, ufl[units,cc], ufl2[units,cc]);
  str(df);
}
