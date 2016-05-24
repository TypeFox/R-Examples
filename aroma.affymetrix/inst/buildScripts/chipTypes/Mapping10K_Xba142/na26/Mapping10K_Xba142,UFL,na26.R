if (interactive()) savehistory();
library("aroma.affymetrix");
log <- Verbose(threshold=-50, timestamp=TRUE);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Settings
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
naVersion <- "26";
user <- "HB";
datestamp <- "20080916";

chipType <- "Mapping10K_Xba142";

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
##  Name: Mapping10K_Xba142
##  Tags: na26,HB20080916
##  Pathname: annotationData/chipTypes/Mapping10K_Xba142/Mapping10K_Xba142,na26,HB20080916.ufl
##  File size: 20.54kB
##  RAM: 0.00MB
##  Number of data rows: 10208
##  File format: v1
##  Dimensions: 10208x1
##  Column classes: integer
##  Number of bytes per column: 2
##  Footer: <platform>Affymetrix</platform><chipType>Mapping10K_Xba142</chipType><createdOn>20080916 18:49:22 PDT</createdOn><createdBy><fullname>Henrik Bengtsson</fullname><email>[...]</email></createdBy><srcFiles><srcFile1><filename>Mapping10K_Xba142.cdf</filename><filesize>9995717</filesize><checksum>ad6ef2568ad7c629c4218b3f402a5cf4</checksum></srcFile1><srcFile2><filename>Mapping10K_Xba142.na26.annot.csv</filename><filesize>16140044</filesize><checksum>666197b57c567684d370cfa8e68057d2</checksum></srcFile2></srcFiles>
##  Chip type: Mapping10K_Xba142
##  Platform: Affymetrix

##                 snp cnp affxSnp other total
##  enzyme1-only 10094   0       0     0 10094
##  missing        110   0       0     4   114
##  total        10204   0       0     4 10208
