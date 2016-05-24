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

chipType <- "Mapping10K_Xba142";

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
##  Name: Mapping10K_Xba142
##  Tags: na26,HB20080916
##  Pathname: annotationData/chipTypes/Mapping10K_Xba142/Mapping10K_Xba142,na26,HB20080916.ugp
##  File size: 50.45kB
##  RAM: 0.00MB
##  Number of data rows: 10208
##  File format: v1
##  Dimensions: 10208x2
##  Column classes: integer, integer
##  Number of bytes per column: 1, 4
##  Footer: <createdOn>20080916 18:46:43 PDT</createdOn><platform>Affymetrix</platform><chipType>Mapping10K_Xba142</chipType><createdBy><fullname>Henrik Bengtsson</fullname><email>[...]</email></createdBy><srcFiles><srcFile1><filename>Mapping10K_Xba142.cdf</filename><filesize>9995717</filesize><checksum>ad6ef2568ad7c629c4218b3f402a5cf4</checksum></srcFile1><srcFile2><filename>Mapping10K_Xba142.na26.annot.csv</filename><filesize>16140044</filesize><checksum>666197b57c567684d370cfa8e68057d2</checksum></srcFile2></srcFiles>
##  Chip type: Mapping10K_Xba142
##  Platform: Affymetrix
##   chromosome        position
##   Min.   :  1.000   Min.   :    93683
##   1st Qu.:  4.000   1st Qu.: 33818850
##   Median :  8.000   Median : 70698761
##   Mean   :  8.811   Mean   : 79654836
##   3rd Qu.: 13.000   3rd Qu.:115978598
##   Max.   : 23.000   Max.   :246860369
##   NA's   :113.000   NA's   :      113

##    1   2   3   4   5   6   7   8   9  10  11  12  13  14  15
##  781 843 708 696 702 694 514 473 487 532 551 480 438 357 284
##   16  17  18  19  20  21  22  23
##  226 167 306 149 204 164  76 263