if (interactive()) savehistory();
library("aroma.affymetrix");
log <- Verbose(threshold=-10, timestamp=TRUE);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Settings
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
naVersion <- "26";
user <- "HB";
datestamp <- "20080915";

chipType <- "Mapping250K_Nsp";

footer <- list(
  createdOn = format(Sys.time(), "%Y%m%d %H:%M:%S", usetz=TRUE),
  createdBy = list(
    fullname = "Henrik Bengtsson", 
    email = sprintf("%s@%s", "henrik.bengtsson", "aroma-project.org")
  ),
  srcFiles = list()
);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Setup required annotation files
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
if (!exists("cdf")) {
  cdf <- AffymetrixCdfFile$byChipType(chipType);
  rm(csv);
}
print(cdf);

if (!exists("csv")) {
  tags <- sprintf(".na%s", naVersion);
  pathname <- AffymetrixNetAffxCsvFile$findByChipType(chipType, tags=tags);
  if (isFile(pathname)) {
    csv <- AffymetrixNetAffxCsvFile(pathname);
  }
  rm(tags);
}
print(csv);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Import UGP from CSV files
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
tags <- sprintf("na%s,%s%s", naVersion, user, datestamp);
ugp <- NULL;
tryCatch({
  ugp <- AromaUgpFile$byChipType(getChipType(cdf), tags=tags);
}, error = function(ex) {})
if (is.null(ugp)) {
  ugp <- AromaUgpFile$allocateFromCdf(cdf, tags=tags);
}
print(ugp);


units <- importFrom(ugp, csv, verbose=log);
str(units);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Update the file footer
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
if (!exists("srcFileTags", mode="list")) {
  srcFileTags <- list();
  srcFiles <- list(cdf, csv);
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

footer <- readFooter(ugp);
footer$createdOn <- format(Sys.time(), "%Y%m%d %H:%M:%S", usetz=TRUE);
footer$createdBy = list(
  fullname = "Henrik Bengtsson", 
  email = sprintf("%s@%s", "henrik.bengtsson", "aroma-project.org")
);
names(srcFileTags) <- sprintf("srcFile%d", seq(along=srcFileTags));
footer$srcFiles <- srcFileTags;
writeFooter(ugp, footer);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Statistics
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
print(ugp);
## AromaUgpFile:
## Name: Mapping250K_Nsp
## Tags: na26,HB20080915
## Pathname: annotationData/chipTypes/Mapping250K_Nsp/Mapping250K_Nsp,na26,HB20080915.ugp
## File size: 1.25MB
## RAM: 0.00MB
## Number of data rows: 262338
## File format: v1
## Dimensions: 262338x2
## Column classes: integer, integer
## Number of bytes per column: 1, 4
## Footer: <createdOn>20080915 22:53:15 PDT</createdOn><platform>Affymetrix</platform><chipType>Mapping250K_Nsp</chipType><createdBy><fullname>Henrik Bengtsson</fullname><email>[...]</email></createdBy><srcFiles><srcFile1><filename>Mapping250K_Nsp.cdf</filename><filesize>194455495</filesize><checksum>59ae263311a2cf63b8d1b9b4cc7d663b</checksum></srcFile1><srcFile2><filename>Mapping250K_Nsp.na26.annot.csv</filename><filesize>447194757</filesize><checksum>5f651b351a0d97e5e9e657f123e6ed05</checksum></srcFile2></srcFiles>
## Chip type: Mapping250K_Nsp
## Platform: Affymetrix

print(summary(ugp));
##  chromosome        position
##  Min.   :  1.000   Min.   :    17408
##  1st Qu.:  4.000   1st Qu.: 32574797
##  Median :  8.000   Median : 70596240
##  Mean   :  8.758   Mean   : 79224244
##  3rd Qu.: 13.000   3rd Qu.:114776301
##  Max.   : 23.000   Max.   :247110269
##  NA's   :775.000   NA's   :      775

print(table(ugp[,1]))
##     1     2     3     4     5     6     7     8     9    10
## 19810 22178 18347 19016 17133 17097 13912 14820 11899 14241
##    11    12    13    14    15    16    17    18    19    20
## 13254 13026 11094  8165  6982  7005  4830  8136  2661  5823
##    21    22    23
##  3927  2511  5696
