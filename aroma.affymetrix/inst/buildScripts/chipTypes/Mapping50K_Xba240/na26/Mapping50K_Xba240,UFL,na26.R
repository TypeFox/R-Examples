if (interactive()) savehistory();
library("aroma.affymetrix");
log <- Verbose(threshold=-50, timestamp=TRUE);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Settings
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
naVersion <- "26";
user <- "HB";
datestamp <- "20080916";

chipType <- "Mapping50K_Xba240";

enzyme <- c("Mapping50K_Hind240"="HindIII", "Mapping50K_Xba240"="XbaI")[chipType];
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
##  Name: Mapping50K_Xba240
##  Tags: na26,HB20080916
##  Pathname: annotationData/chipTypes/Mapping50K_Xba240/Mapping50K_Xba240,na26,HB20080916.ufl
##  File size: 115.87kB
##  RAM: 0.00MB
##  Number of data rows: 59015
##  File format: v1
##  Dimensions: 59015x1
##  Column classes: integer
##  Number of bytes per column: 2
##  Footer: <platform>Affymetrix</platform><chipType>Mapping50K_Xba240</chipType><createdOn>20080916 18:35:23 PDT</createdOn><createdBy><fullname>Henrik Bengtsson</fullname><email>[...]</email></createdBy><srcFiles><srcFile1><filename>Mapping50K_Xba240.CDF</filename><filesize>57703961</filesize><checksum>7f026f6472f7721255717fb3f453be61</checksum></srcFile1><srcFile2><filename>Mapping50K_Xba240.na26.annot.csv</filename><filesize>94903244</filesize><checksum>eab081f8e918293016a760c5714d1198</checksum></srcFile2></srcFiles>
##  Chip type: Mapping50K_Xba240
##  Platform: Affymetrix

##                 snp cnp affxSnp other total
##  enzyme1-only 58616   0       0     0 58616
##  missing        344   0       0    55   399
##  total        58960   0       0    55 59015
