if (interactive()) savehistory();
library("aroma.affymetrix");
log <- Verbose(threshold=-10, timestamp=TRUE);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Settings
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
naVersion <- "26";
user <- "HB";
datestamp <- "20080821";

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
## Tags: Full,na26,HB20080821
## Pathname: annotationData/chipTypes/GenomeWideSNP_6/
##           GenomeWideSNP_6,Full,na26,HB20080821.ugp
## File size: 8.97MB
## RAM: 0.00MB
## Number of data rows: 1881415
## File format: v1
## Dimensions: 1881415x2
## Column classes: integer, integer
## Number of bytes per column: 1, 4
## Footer: ...
## Chip type: GenomeWideSNP_6,Full
## Platform: Affymetrix
getChromosomeStats(ugp);
## GenomeWideSNP_6,Full,na26,HB20080821:


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# WHAT'S NEW:
#
# o na24 -> na26
#   Only minor modifications for non-missing values:
#   - three loci changed chromosomes
#   - an additional 23 loci changed positions, of which only 17 moved
#     more than 2 base pairs.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
ugp <- AromaUgpFile$byChipType("GenomeWideSNP_6,Full", tags="na26");
ugp2 <- AromaUgpFile$byChipType("GenomeWideSNP_6,Full", tags="na24");
cdf <- AffymetrixCdfFile$byChipType(getChipType(ugp));

print(table(ugp[,1], exclude=NULL));
##       1      2      3      4      5      6      7      8
##  146401 153663 127766 120296 115672 112825 100996  98277
##       9     10     11     12     13     14     15     16
##   82168  93592  89525  87321  66067  57103  53556  54182
##      17     18     19     20     21     22     23     24
##   46632  52093  30299  43628  25111  24484  87198   9485
##      25   <NA>
##     445   2630

print(table(ugp2[,1], exclude=NULL));
##       1      2      3      4      5      6      7      8
##  146524 153732 127815 120360 115731 112895 101093  98306
##       9     10     11     12     13     14     15     16
##   82225  93655  89615  87372  66106  57121  53595  54215
##      17     18     19     20     21     22     23     24
##   46678  52109  30362  43648  25129  24513  87204   9486
##    <NA>
##    1926

## Differences on ChrX and ChrY
for (cc in 23:24) {
  units <- getUnitsAt(ugp, cc);
  units2 <- getUnitsAt(ugp2, cc);
  du <- sort(unique(c(setdiff(units, units2), setdiff(units2, units))));
  df <- data.frame(unit=du, unitName=getUnitNames(cdf, units=du), inNew=is.element(du, units), inOld=is.element(du, units2));
  print(df);
}

## ChrX:
##       unit      unitName inNew inOld
##  1   62226 SNP_A-4213851 FALSE  TRUE
##  2   71775 SNP_A-4270117 FALSE  TRUE
##  3   90275 SNP_A-2036574 FALSE  TRUE
##  4  107326 SNP_A-2056054 FALSE  TRUE
##  5  150982 SNP_A-4223523  TRUE FALSE
##  6  263684 SNP_A-4295551 FALSE  TRUE
##  7  280319 SNP_A-4297740 FALSE  TRUE
##  8  287813 SNP_A-2261398 FALSE  TRUE
##  9  296412 SNP_A-2271143 FALSE  TRUE
##  10 340487 SNP_A-1786773 FALSE  TRUE
##  11 348472 SNP_A-1795905 FALSE  TRUE
##  12 358905 SNP_A-1807644 FALSE  TRUE
##  13 361602 SNP_A-1810640 FALSE  TRUE
##  14 365893 SNP_A-1815384 FALSE  TRUE
##  15 370741 SNP_A-1820836  TRUE FALSE
##  16 372245 SNP_A-4246284  TRUE FALSE
##  17 423139 SNP_A-1880658 FALSE  TRUE
##  18 457297 SNP_A-1919606 FALSE  TRUE
##  19 459820 SNP_A-4204363 FALSE  TRUE
##  20 480322 SNP_A-1945393 FALSE  TRUE
##  21 494913 SNP_A-4207824 FALSE  TRUE
##  22 538245 SNP_A-8288469 FALSE  TRUE
##  23 566478 SNP_A-8365290 FALSE  TRUE
##  24 576845 SNP_A-8370182 FALSE  TRUE
##  25 581409 SNP_A-8377813  TRUE FALSE
##  26 584896 SNP_A-8378063  TRUE FALSE
##  27 606965 SNP_A-8458950  TRUE FALSE
##  28 608486 SNP_A-8461206  TRUE FALSE
##  29 631792 SNP_A-8537736  TRUE FALSE
##  30 637955 SNP_A-8534336  TRUE FALSE
##  31 658300 SNP_A-8565568  TRUE FALSE
##  32 668862 SNP_A-8588386 FALSE  TRUE
##  33 704564 SNP_A-8631242 FALSE  TRUE
##  34 721717 SNP_A-8651867  TRUE FALSE
##  35 743165 SNP_A-8682714  TRUE FALSE
##  36 744586 SNP_A-8682083  TRUE FALSE
##  37 758776 SNP_A-8702321 FALSE  TRUE
##  38 763415 SNP_A-8696584  TRUE FALSE
##  39 787360 SNP_A-8357497  TRUE FALSE
##  40 823570 SNP_A-8489739  TRUE FALSE
##  41 870048 SNP_A-8422929 FALSE  TRUE
##  42 894880 SNP_A-8492366  TRUE FALSE
##  43 920930 SNP_A-8531829  TRUE FALSE
##  44 932378 SNP_A-8284454  TRUE FALSE

## ChrY:
##       unit      unitName inNew inOld
##  1  191853 SNP_A-2152277 FALSE  TRUE
##  2  495041 SNP_A-8290289  TRUE FALSE
##  3  495193 SNP_A-8477472 FALSE  TRUE
##  4  495229 SNP_A-8539202  TRUE FALSE
##  5  495230 SNP_A-8539203  TRUE FALSE
##  6  495438 SNP_A-8647645 FALSE  TRUE
##  7  495472 SNP_A-8677558 FALSE  TRUE
##  8  495529 SNP_A-8714042 FALSE  TRUE
##  9  495805 SNP_A-8511463  TRUE FALSE
##  10 495806 SNP_A-8511470  TRUE FALSE
##  11 495865 SNP_A-8528164 FALSE  TRUE
##  12 702145 SNP_A-8629588 FALSE  TRUE
##  13 710044 SNP_A-8638351  TRUE FALSE


for (cc in 1:nbrOfColumns(ugp)) {
  units <- whichVector(ugp[,cc] != ugp2[,cc], na.rm=FALSE);
  df <- cbind(units, ugp[units,cc], ugp2[units,cc], ugp[units,cc]-ugp2[units,cc]);
  colnames(df)[ncol(df)] <- "delta";
  print(df);
}

##    chromosome chromosome delta
##  1         15         22    -7
##  2          7         24   -17
##  3          1          3    -2

##      position  position     delta
##  1   20016405  47401631 -27385226
##  2  149942827  57728335  92214492
##  3  112154796 112154795         1
##  4   85836181  82121726   3714455
##  5   25697318  25697317         1
##  6  148400695  78226876  70173819
##  7   72582965  66455575   6127390
##  8  141510003  75974344  65535659
##  9   45204796  28233161  16971635
##  10  88282234  88282233         1
##  11  94370183  94370182         1
##  12  26768168   5462681  21305487
##  13 100234698  82328615  17906083
##  14  20408524  20452780    -44256
##  15 229505938 152760680  76745258
##  16  50070297   8427474  41642823
##  17  70958284  62600568   8357716
##  18  47156058  47156154       -96
##  19  39862439  39862440        -1
##  20  69852830  63620439   6232391
##  21  94700194  78792665  15907529
##  22  29072915  29072913         2
##  23 196317488 161362741  34954747
##  24  53493738  40654350  12839388
##  25  25595288   7749030  17846258
##  26 177588588 177588587         1
