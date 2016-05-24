if (interactive()) savehistory();
library("aroma.affymetrix");
log <- Verbose(threshold=-10, timestamp=TRUE);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Settings
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
naVersion <- "28";
user <- "HB";
datestamp <- "20090519";

chipType <- "Cytogenetics_Array";

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
  rm(csvList);
}
print(cdf);

if (!exists("csvList", mode="list")) {
  csvList <- list();

  tagsList <- c(
      main=sprintf(".na%s", naVersion),
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
  ## Cytogenetics_Array.na28.annot.csv:    int [1:400103] 325823 117191
  ## Cytogenetics_Array.cn.na28.annot.csv: int [1:2387593] 2771313 2771325
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
AromaUgpFile:
Name: Cytogenetics_Array
Tags: na28,HB20090519
Full name: Cytogenetics_Array,na28,HB20090519
Pathname: annotationData/chipTypes/Cytogenetics_Array/Cytogenetics_Array,na28,HB20090519.ugp
File size: 13.33 MB (13979333 bytes)
RAM: 0.00 MB
Number of data rows: 2795708
File format: v1
Dimensions: 2795708x2
Column classes: integer, integer
Number of bytes per column: 1, 4
Footer: <createdOn>20090519 23:15:52 PDT</createdOn><platform>Affymetrix</platform><chipType>Cytogenetics_Array</chipType><createdBy><fullname>Henrik Bengtsson</fullname><email>[...]</email></createdBy><srcFiles><srcFile1><filename>Cytogenetics_Array.CDF</filename><filesize>564465214</filesize><checksum>17ea0ac70c196f20a1825dff8b5e773c</checksum></srcFile1><srcFile2><filename>Cytogenetics_Array.na28.annot.csv</filename><filesize>141110783</filesize><checksum>eeca051fc6c7822173e276cac14e0ad4</checksum></srcFile2><srcFile3><filename>Cytogenetics_Array.cn.na28.annot.csv</filename><filesize>852047697</filesize><checksum>63490620be89e306767de896e883b5fa</checksum></srcFile3></srcFiles>
Chip type: Cytogenetics_Array
Platform: Affymetrix

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# WHAT'S NEW:
#
# o na28
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
ugp <- AromaUgpFile$byChipType("Cytogenetics_Array", tags="na28");

options(width=60);
print(table(ugp[,1], exclude=NULL));
     1      2      3      4      5      6      7      8
222413 226698 192608 169326 167528 166405 155857 139341
     9     10     11     12     13     14     15     16
121665 132176 130925 124991  89622  86105  83585  81726
    17     18     19     20     21     22     23     24
 84215  72954  59305  60120  33911  37351 132830  14386
    25   <NA>
  1653   8012

