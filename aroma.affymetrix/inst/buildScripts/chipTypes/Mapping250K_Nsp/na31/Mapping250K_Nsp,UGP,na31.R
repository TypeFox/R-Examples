if (interactive()) savehistory();
library("aroma.affymetrix");
log <- Verbose(threshold=-10, timestamp=TRUE);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Settings
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
naVersion <- "31";
user <- "HB";
datestamp <- "20101007";

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
## Tags: na31,HB20101007
## Full name: Mapping250K_Nsp,na31,HB20101007
## Pathname: annotationData/chipTypes/Mapping250K_Nsp/Mapping250K_Nsp,na31,HB20101007.ugp
## File size: 1.25 MB (1312308 bytes)
## RAM: 0.01 MB
## Number of data rows: 262338
## File format: v1
## Dimensions: 262338x2
## Column classes: integer, integer
## Number of bytes per column: 1, 4
## Footer: <createdOn>20101007 14:59:46 PDT</createdOn><platform>Affymetrix</platform><chipType>Mapping250K_Nsp</chipType><createdBy><fullname>Henrik Bengtsson</fullname><email>[...]</email></createdBy><srcFiles><srcFile1><filename>Mapping250K_Nsp.cdf</filename><filesize>194455495</filesize><checksum>59ae263311a2cf63b8d1b9b4cc7d663b</checksum></srcFile1><srcFile2><filename>Mapping250K_Nsp.na31.annot.csv</filename><filesize>400278099</filesize><checksum>dbd9299c348e1b688f752343dc9e8e57</checksum></srcFile2></srcFiles>
## Chip type: Mapping250K_Nsp
## Platform: Affymetrix

print(summary(ugp));
##  chromosome        position
##  Min.   :  1.000   Min.   :    18901
##  1st Qu.:  4.000   1st Qu.: 32979020
##  Median :  8.000   Median : 71259947
##  Mean   :  8.759   Mean   : 79599960
##  3rd Qu.: 13.000   3rd Qu.:114973065
##  Max.   : 24.000   Max.   :249143646
##  NA's   :357.000   NA's   :      357

print(table(ugp[,1]))
##     1     2     3     4     5     6     7     8     9    10
## 19858 22206 18359 19053 17153 17106 13932 14835 11942 14264
##    11    12    13    14    15    16    17    18    19    20
## 13290 13047 11075  8175  7003  7010  4844  8145  2691  5839
##    21    22    23    24
##  3936  2499  5714     5
