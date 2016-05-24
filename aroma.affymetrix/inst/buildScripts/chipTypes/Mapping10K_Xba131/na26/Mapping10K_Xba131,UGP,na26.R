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

chipType <- "Mapping10K_Xba131";

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

##  AromaUgpFile:
##  Name: Mapping10K_Xba131
##  Tags: na26,HB20080916
##  Pathname: annotationData/chipTypes/Mapping10K_Xba131/Mapping10K_Xba131,na26,HB20080916.ugp
##  File size: 57.07kB
##  RAM: 0.00MB
##  Number of data rows: 11564
##  File format: v1
##  Dimensions: 11564x2
##  Column classes: integer, integer
##  Number of bytes per column: 1, 4
##  Footer: <createdOn>20080916 18:51:07 PDT</createdOn><platform>Affymetrix</platform><chipType>Mapping10K_Xba131</chipType><createdBy><fullname>Henrik Bengtsson</fullname><email>[...]</email></createdBy><srcFiles><srcFile1><filename>Mapping10K_Xba131.cdf</filename><filesize>11311508</filesize><checksum>a18490aadbdc6407332ac3bb12d55a75</checksum></srcFile1><srcFile2><filename>Mapping10K_Xba131.na26.annot.csv</filename><filesize>18345600</filesize><checksum>4d21244a2f48e0484c01cbf50659a9ef</checksum></srcFile2></srcFiles>
##  Chip type: Mapping10K_Xba131
##  Platform: Affymetrix
##   chromosome        position
##   Min.   :  1.000   Min.   :    93683
##   1st Qu.:  4.000   1st Qu.: 33855043
##   Median :  8.000   Median : 70699090
##   Mean   :  8.742   Mean   : 79648642
##   3rd Qu.: 13.000   3rd Qu.:115969408
##   Max.   : 23.000   Max.   :246860369
##   NA's   :132.000   NA's   :      132

##    1   2   3   4   5   6   7   8   9  10  11  12  13  14  15
##  890 960 814 808 779 791 584 555 541 610 638 545 490 400 329
##   16  17  18  19  20  21  22  23
##  260 188 345  97 223 195  81 309