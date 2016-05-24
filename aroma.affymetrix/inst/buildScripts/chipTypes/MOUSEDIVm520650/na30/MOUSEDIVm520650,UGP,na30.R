if (interactive()) savehistory();
library("aroma.affymetrix");
verbose <- Verbose(threshold=-10, timestamp=TRUE);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Settings
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
naVersion <- "30";
genomeVersion <- "mm9";
user <- "HB";
datestamp <- "20100603";

chipType <- "MOUSEDIVm520650";


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Setup required annotation files
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
if (!exists("cdf")) {
  cdf <- AffymetrixCdfFile$byChipType(chipType);
  rm(csvList);
}
print(cdf);

if (!exists("csvList", mode="list")) {
  csvList <- list();

  tagsList <- c(
      main=sprintf(".na%s", naVersion),
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
# Setup/Allocate UGP
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


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Import from CSV files
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
for (kk in seq(along=csvList)) {
  csv <- csvList[[kk]];
  print(csv);
  units <- importFrom(ugp, csv, verbose=verbose);
  str(units);
  ## *.na30.annot.csv:  int [1:623124] 800990 808457 907857 ...
  ## *.cn.na30.annot.csv: int [1:1832538] 1542481 1542480 1542479 ...
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
footer$createdBy <- list(
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
## Name: MOUSEDIVm520650
## Tags: na30,mm9,HB20190603
## Full name: MOUSEDIVm520650,na30,mm9,HB20190603
## Pathname: annotationData/chipTypes/MOUSEDIVm520650/MOUSEDIVm520650,na30,mm9,HB20190603.ugp
## File size: 11.72 MB (12294266 bytes)
## RAM: 0.01 MB
## Number of data rows: 2458697
## File format: v1
## Dimensions: 2458697x2
## Column classes: integer, integer
## Number of bytes per column: 1, 4
## Footer: <createdOn>20100603 05:49:21 PDT</createdOn><platform>Affymetrix</platform><chipType>MOUSEDIVm520650</chipType><createdBy><fullname>Henrik Bengtsson</fullname><email>[...]</email></createdBy><srcFiles><srcFile1><filename>MOUSEDIVm520650.CDF</filename><filesize>668362285</filesize><checksum>1d83d7dde5e4816f4be315f6863f6a7c</checksum></srcFile1><srcFile2><filename>MOUSEDIVm520650.na30.annot.csv</filename><filesize>358437027</filesize><checksum>e3ef35245954d5a08ae299b95666fc3a</checksum></srcFile2><srcFile3><filename>MOUSEDIVm520650.cn.na30.annot.csv</filename><filesize>795251151</filesize><checksum>5b5bbdb60ed443b2616358beea06ddaf</checksum></srcFile3></srcFiles>
## Chip type: MOUSEDIVm520650
## Platform: Affymetrix


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# WHAT'S NEW:
#
# o ... -> na30
#   - First one for this chip type.  We note that SNP CSV file uses "M"
#     for mitochondrial DNA, whereas the CN CSV uses "MT".  This import
#     script recognizes both.  The CN CSV also contains 227,740 "NA":s,
#     which is interpreted as "unknown" locations.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
ugp <- AromaUgpFile$byChipType("MOUSEDIVm520650", tags="na30");

options(width=60);
print(table(ugp[,1], exclude=NULL));

##      1      2      3      4      5      6      7      8
## 136405 170236 107112 142492 133169 108251 171812 100077
##      9     10     11     12     13     14     15     16
## 106036  88738 135094 108539  89828 108188  79179  64466
##     17     18     19     23     24     25   <NA>
##  90052  58269  63067 165816    987    109 230775
