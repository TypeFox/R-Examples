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


for (kk in seq(along=csvList)) {
  csv <- csvList[[kk]];
  print(csv);
  units <- importFrom(ugp, csv, verbose=log);
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
## Name: GenomeWideSNP_6
## Tags: Full,na27.1,HB20090519
## Full name: GenomeWideSNP_6,Full,na27.1,HB20090519
## Pathname: annotationData/chipTypes/GenomeWideSNP_6/GenomeWideSNP_6,Full,na27.1,HB20090519.ugp
## File size: 8.97 MB (9407869 bytes)
## RAM: 0.01 MB
## Number of data rows: 1881415
## File format: v1
## Dimensions: 1881415x2
## Column classes: integer, integer
## Number of bytes per column: 1, 4
## Footer: <createdOn>20090519 18:56:13 PDT</createdOn><platform>Affymetrix</platform><chipType>GenomeWideSNP_6,Full</chipType><createdBy><fullname>Henrik Bengtsson</fullname><email>[...]</email></createdBy><srcFiles><srcFile1><filename>GenomeWideSNP_6,Full.CDF</filename><filesize>493291745</filesize><checksum>3fbe0f6e7c8a346105238a3f3d10d4ec</checksum></srcFile1><srcFile2><filename>GenomeWideSNP_6.na27.1.annot.csv</filename><filesize>1473437359</filesize><checksum>2b1cc87850a4c5e762b52c157689f570</checksum></srcFile2><srcFile3><filename>GenomeWideSNP_6.cn.na27.annot.csv</filename><filesize>490970502</filesize><checksum>2ad99a80d26b0af8f14b1a5e6d92d81c</checksum></srcFile3></srcFiles>
## Chip type: GenomeWideSNP_6,Full
## Platform: Affymetrix
getChromosomeStats(ugp);
## GenomeWideSNP_6,Full,na26,HB20080821:


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# WHAT'S NEW:
#
# o na26 -> na27
#   Two units (932039, 1872834) where moved from ChrX to ChrY.
#   Same location.
# o na24 -> na26
#   Only minor modifications for non-missing values:
#   - three loci changed chromosomes
#   - an additional 23 loci changed positions, of which only 17 moved
#     more than 2 base pairs.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
ugp <- AromaUgpFile$byChipType("GenomeWideSNP_6,Full", tags="na27.1");
ugp2 <- AromaUgpFile$byChipType("GenomeWideSNP_6,Full", tags="na26");
cdf <- AffymetrixCdfFile$byChipType(getChipType(ugp));

options(width=60);
print(table(ugp[,1], exclude=NULL));
##       1      2      3      4      5      6      7      8
##  146401 153663 127766 120296 115672 112825 100996  98277
##       9     10     11     12     13     14     15     16
##   82168  93592  89525  87321  66067  57103  53556  54182
##      17     18     19     20     21     22     23     24
##   46632  52093  30299  43628  25111  24484  87200   9483
##      25   <NA>
##     445   2630

print(table(ugp2[,1], exclude=NULL));
##       1      2      3      4      5      6      7      8
##  146401 153663 127766 120296 115672 112825 100996  98277
##       9     10     11     12     13     14     15     16
##   82168  93592  89525  87321  66067  57103  53556  54182
##      17     18     19     20     21     22     23     24
##   46632  52093  30299  43628  25111  24484  87198   9485
##      25   <NA>
##     445   2630

rr <- whichVector(ugp[,1] != ugp2[,1]);
print(rr);
## [1] 932039 1872834;
print(cbind(ugp[rr,], ugp2[rr,]));
##   chromosome position chromosome position
## 1         23   169542         24   169542
## 2         23  2675210         24  2675210

rr <- whichVector(ugp[,2] != ugp2[,2]);
print(rr);
## integer(0)
