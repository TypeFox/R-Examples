if (interactive()) savehistory();
library("aroma.affymetrix");
log <- Verbose(threshold=-50, timestamp=TRUE);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Settings
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
naVersion <- "26";
user <- "HB";
datestamp <- "20080916";

chipType <- "Mapping10K_Xba131";

enzyme <- "XbaI";
print(enzyme);

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
# Import UFL from CSV files
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
tags <- sprintf("na%s,%s%s", naVersion, user, datestamp);
ufl <- NULL;
tryCatch({
  ufl <- AromaUflFile$byChipType(getChipType(cdf), tags=tags);
}, error = function(ex) {})
if (is.null(ufl)) {
  ufl <- AromaUflFile$allocateFromCdf(cdf, tags=tags);
}
print(ufl);


stopifnot(!is.na(enzyme));
units <- importFrom(ufl, csv, enzymes=enzyme, verbose=log);
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

footer <- readFooter(ufl);
footer$createdOn <- format(Sys.time(), "%Y%m%d %H:%M:%S", usetz=TRUE);
footer$createdBy = list(
  fullname = "Henrik Bengtsson", 
  email = sprintf("%s@%s", "henrik.bengtsson", "aroma-project.org")
);
names(srcFileTags) <- sprintf("srcFile%d", seq(along=srcFileTags));
footer$srcFiles <- srcFileTags;
writeFooter(ufl, footer);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Statistics
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
print(ufl);
x <- summaryOfUnits(ufl);
print(x);

##  AromaUflFile:
##  Name: Mapping10K_Xba131
##  Tags: na26,HB20080916
##  Pathname: annotationData/chipTypes/Mapping10K_Xba131/Mapping10K_Xba131,na26,HB20080916.ufl
##  File size: 23.19kB
##  RAM: 0.00MB
##  Number of data rows: 11564
##  File format: v1
##  Dimensions: 11564x1
##  Column classes: integer
##  Number of bytes per column: 2
##  Footer: <platform>Affymetrix</platform><chipType>Mapping10K_Xba131</chipType><createdOn>20080916 18:52:16 PDT</createdOn><createdBy><fullname>Henrik Bengtsson</fullname><email>[...]</email></createdBy><srcFiles><srcFile1><filename>Mapping10K_Xba131.cdf</filename><filesize>11311508</filesize><checksum>a18490aadbdc6407332ac3bb12d55a75</checksum></srcFile1><srcFile2><filename>Mapping10K_Xba131.na26.annot.csv</filename><filesize>18345600</filesize><checksum>4d21244a2f48e0484c01cbf50659a9ef</checksum></srcFile2></srcFiles>
##  Chip type: Mapping10K_Xba131
##  Platform: Affymetrix

##                 snp cnp affxSnp other total
##  enzyme1-only 11431   0       0     0 11431
##  missing        129   0       0     4   133
##  total        11560   0       0     4 11564
