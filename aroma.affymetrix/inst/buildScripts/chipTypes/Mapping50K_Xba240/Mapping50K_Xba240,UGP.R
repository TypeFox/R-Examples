if (interactive()) savehistory();
library("aroma.affymetrix");
library("R.menu");

verbose <- Verbose(threshold=-10, timestamp=TRUE);
options(width=60);

chipType <- "Mapping50K_Xba240";


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# User settings
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
## setOption(aromaSettings, "user/initials", "HB");
## setOption(aromaSettings, "user/fullname", "Henrik Bengtsson");
## email <- sprintf("%s@%s", "henrik.bengtsson", "aroma-project.org")
## setOption(aromaSettings, "user/email", email);
## saveAnywhere(aromaSettings);

fullname <- getOption(aromaSettings, "user/fullname");
stopifnot(!is.null(fullname));
email <- getOption(aromaSettings, "user/email");
stopifnot(!is.null(email));
user <- getOption(aromaSettings, "user/initials");
stopifnot(!is.null(user));


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Settings
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
genomeVersions <- c("30"="hg18", "31"="hg19");
naVersions <- names(genomeVersions);
choices <- sprintf("na%s (%s)", naVersions, genomeVersions);
choice <- textMenu(choices, title="Choose NetAffx version: ", value=FALSE);
naVersion <- naVersions[choice];
genomeVersion <- genomeVersions[naVersion];
datestamp <- format(Sys.Date(), format="%Y%m%d");



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Setup required annotation files
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
cdf <- AffymetrixCdfFile$byChipType(chipType);
print(cdf);

tags <- sprintf(".na%s", naVersion);
pathname <- AffymetrixNetAffxCsvFile$findByChipType(chipType, tags=tags);
if (isFile(pathname)) {
  csv <- AffymetrixNetAffxCsvFile(pathname);
}
rm(tags);
print(csv);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Import UGP from CSV files
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
tags <- sprintf("na%s,%s,%s%s", naVersion, genomeVersion, user, datestamp);
ugp <- NULL;
tryCatch({
  ugp <- AromaUgpFile$byChipType(getChipType(cdf), tags=tags);
}, error = function(ex) {})
if (is.null(ugp)) {
  ugp <- AromaUgpFile$allocateFromCdf(cdf, tags=tags);
}
print(ugp);


units <- importFrom(ugp, csv, verbose=verbose);
str(units);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Update the file footer
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
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

footer <- readFooter(ugp);
footer$createdOn <- format(Sys.time(), "%Y%m%d %H:%M:%S", usetz=TRUE);
footer$createdBy = list(
  fullname = fullname, 
  email = email
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
## Name: Mapping50K_Xba240
## Tags: na31,hg19,HB20110328
## Full name: Mapping50K_Xba240,na31,hg19,HB20110328
## Pathname: annotationData/chipTypes/Mapping50K_Xba240/Mapping50K_Xba240,na31,hg19
## ,HB20110328.ugp
## File size: 288.77 kB (295699 bytes)
## RAM: 0.01 MB
## Number of data rows: 59015
## File format: v1
## Dimensions: 59015x2
## Column classes: integer, integer
## Number of bytes per column: 1, 4
## Footer: <createdOn>20110328 18:46:23 PDT</createdOn><platform>Affymetrix</platfo
## rm><chipType>Mapping50K_Xba240</chipType><createdBy><fullname>Henrik Bengtsson</
## fullname><email>[...]</email></createdBy><srcFiles><srcFile1><fil
## ename>Mapping50K_Xba240.CDF</filename><filesize>57703961</filesize><checksum>7f0
## 26f6472f7721255717fb3f453be61</checksum></srcFile1><srcFile2><filename>Mapping50
## K_Xba240.na31.annot.csv</filename><filesize>84705402</filesize><checksum>4cc4997
## 35dd13177331994cde4e9f323</checksum></srcFile2></srcFiles>
## Chip type: Mapping50K_Xba240
## Platform: Affymetrix
## 
##  chromosome        position
##  Min.   :  1.000   Min.   :   110819
##  1st Qu.:  4.000   1st Qu.: 35046332
##  Median :  7.000   Median : 72991749
##  Mean   :  8.504   Mean   : 80414591
##  3rd Qu.: 12.000   3rd Qu.:114933994
##  Max.   : 24.000   Max.   :248818466
##  NA's   :485.000   NA's   :      485
## 
##    1    2    3    4    5    6    7    8    9   10   11   12
## 4662 5271 3863 4227 4141 4094 3609 3422 2436 2948 2892 2684
##   13   14   15   16   17   18   19   20   21   22   23   24
## 2555 2081 1587 1238  968 1835  363 1100 1023  326 1203    2
