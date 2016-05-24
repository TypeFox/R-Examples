if (interactive()) savehistory();
library("aroma.affymetrix");
log <- Verbose(threshold=-10, timestamp=TRUE);
options(width=60);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Settings
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
naVersion <- "26";
user <- "HB";
datestamp <- "20080916";

chipType <- "Mapping50K_Hind240";

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
print(summary(ugp));
print(table(ugp[,1]));

## AromaUgpFile:
## Name: Mapping50K_Hind240
## Tags: na26,HB20080916
## Pathname: annotationData/chipTypes/Mapping50K_Hind240/Mapping50K_Hind240,na26,HB20080916.ugp
## File size: 280.39kB
## RAM: 0.00MB
## Number of data rows: 57299
## File format: v1
## Dimensions: 57299x2
## Column classes: integer, integer
## Number of bytes per column: 1, 4
## Footer: <createdOn>20080916 18:11:33 PDT</createdOn><platform>Affymetrix</platform><chipType>Mapping50K_Hind240</chipType><createdBy><fullname>Henrik Bengtsson</fullname><email>[...]</email></createdBy><srcFiles><srcFile1><filename>Mapping50K_Hind240.CDF</filename><filesize>56029390</filesize><checksum>84239fce7a443dd528e060731a8df77e</checksum></srcFile1><srcFile2><filename>Mapping50K_Hind240.na26.annot.csv</filename><filesize>92433908</filesize><checksum>fc9a48388f38491b5be7ffaf1fdc3861</checksum></srcFile2></srcFiles>
## Chip type: Mapping50K_Hind240
## Platform: Affymetrix

##  chromosome        position
##  Min.   :  1.000   Min.   :    48603
##  1st Qu.:  4.000   1st Qu.: 34667112
##  Median :  7.000   Median : 72677621
##  Mean   :  8.402   Mean   : 80405004
##  3rd Qu.: 12.000   3rd Qu.:114826216
##  Max.   : 23.000   Max.   :246727435
##  NA's   :363.000   NA's   :      363

##    1    2    3    4    5    6    7    8    9   10   11   12
## 4541 5072 3962 4342 4215 3968 3444 3549 2357 2743 2466 2592
##   13   14   15   16   17   18   19   20   21   22   23
## 2661 1931 1440 1145  985 1731  326  993  883  433 1157
