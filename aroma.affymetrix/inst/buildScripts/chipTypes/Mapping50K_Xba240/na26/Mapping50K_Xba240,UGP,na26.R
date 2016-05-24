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

chipType <- "Mapping50K_Xba240";

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
##  Name: Mapping50K_Xba240
##  Tags: na26,HB20080916
##  Pathname: C:/Users/hb/Documents/My Data/annotationData/chipTypes/Mapping50K_Xba2
##  40/Mapping50K_Xba240,na26,HB20080916.ugp
##  File size: 288.77kB
##  RAM: 0.00MB
##  Number of data rows: 59015
##  File format: v1
##  Dimensions: 59015x2
##  Column classes: integer, integer
##  Number of bytes per column: 1, 4
##  Footer: <createdOn>20080916 18:16:09 PDT</createdOn><platform>Affymetrix</platform><chipType>Mapping50K_Xba240</chipType><createdBy><fullname>Henrik Bengtsson</fullname><email>[...]</email></createdBy><srcFiles><srcFile1><filename>Mapping50K_Xba240.CDF</filename><filesize>57703961</filesize><checksum>7f026f6472f7721255717fb3f453be61</checksum></srcFile1><srcFile2><filename>Mapping50K_Xba240.na26.annot.csv</filename><filesize>94903244</filesize><checksum>eab081f8e918293016a760c5714d1198</checksum></srcFile2></srcFiles>
##  Chip type: Mapping50K_Xba240
##  Platform: Affymetrix

##   chromosome        position
##   Min.   :  1.000   Min.   :    93683
##   1st Qu.:  4.000   1st Qu.: 34636629
##   Median :  7.000   Median : 72249739
##   Mean   :  8.507   Mean   : 80010574
##   3rd Qu.: 12.000   3rd Qu.:114666170
##   Max.   : 24.000   Max.   :246885089
##   NA's   :390.000   NA's   :      390

##     1    2    3    4    5    6    7    8    9   10   11   12
##  4669 5274 3864 4231 4149 4111 3612 3422 2439 2953 2896 2685
##    13   14   15   16   17   18   19   20   21   22   23   24
##  2567 2083 1590 1245  971 1837  364 1101 1027  330 1204    1
