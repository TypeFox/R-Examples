if (interactive()) savehistory();
library("aroma.affymetrix");
log <- Verbose(threshold=-10, timestamp=TRUE);
options(width=60);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Settings
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
naVersion <- "26";
user <- "HB";
datestamp <- "20080915";

chipType <- "Mapping250K_Sty";

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
##  AromaUgpFile:
##  Name: Mapping250K_Sty
##  Tags: na26,HB20080915
##  Pathname: annotationData/chipTypes/Mapping250K_Sty/Mapping250K_Sty,na26,HB20080915.ugp
##  File size: 1.14MB
##  RAM: 0.00MB
##  Number of data rows: 238378
##  File format: v1
##  Dimensions: 238378x2
##  Column classes: integer, integer
##  Number of bytes per column: 1, 4
##  Footer: <createdOn>20080915 21:55:14 PDT</createdOn><platform>Affymetrix</platform><chipType>Mapping250K_Sty</chipType><createdBy><fullname>Henrik Bengtsson</fullname><email>[...]</email></createdBy><srcFiles><srcFile1><filename>Mapping250K_Sty.cdf</filename><filesize>184371335</filesize><checksum>8c7ffa88bb444d953214f601bb998e51</checksum></srcFile1><srcFile2><filename>Mapping250K_Sty.na26.annot.csv</filename><filesize>412941971</filesize><checksum>eb1ced7bd45aa7c6f6edf9fb708c7441</checksum></srcFile2></srcFiles>
##  Chip type: Mapping250K_Sty
##  Platform: Affymetrix

print(summary(ugp));
##  chromosome        position
##  Min.   :  1.000   Min.   :     2994
##  1st Qu.:  4.000   1st Qu.: 31306881
##  Median :  8.000   Median : 67082398
##  Mean   :  9.117   Mean   : 77333484
##  3rd Qu.: 13.000   3rd Qu.:114799352
##  Max.   : 23.000   Max.   :247135059
##  NA's   :677.000   NA's   :      677

print(table(ugp[,1]));
##      1     2     3     4     5     6     7     8     9    10
##  20273 19147 15396 13236 14854 14296 11822 12601 10900 14190
##     11    12    13    14    15    16    17    18    19    20
##  12928 11863  8052  7531  7323  8272  6403  6726  3670  6555
##     21    22    23
##   3182  3669  4812
