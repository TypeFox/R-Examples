if (interactive()) savehistory();
library("aroma.affymetrix");
log <- Verbose(threshold=-50, timestamp=TRUE);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Settings
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
naVersion <- "31";
user <- "HB";
datestamp <- "20101007";

chipType <- "Mapping250K_Nsp";

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


enzyme <- gsub("Mapping250K_(Nsp|Sty)", "\\1I", getChipType(ufl));
print(enzyme);
stopifnot(nchar(enzyme) == 4);
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
##  Name: Mapping250K_Nsp
##  Tags: na26,HB20080916
##  Pathname: annotationData/chipTypes/Mapping250K_Nsp/Mapping250K_Nsp,na26,HB20080916.ufl
##  File size: 512.98kB
##  RAM: 0.00MB
##  Number of data rows: 262338
##  File format: v1
##  Dimensions: 262338x1
##  Column classes: integer
##  Number of bytes per column: 2
##  Footer: <platform>Affymetrix</platform><chipType>Mapping250K_Nsp</chipType><createdOn>20080916 17:23:12 PDT</createdOn><createdBy><fullname>Henrik Bengtsson</fullname><email>[...]</email></createdBy><srcFiles><srcFile1><filename>Mapping250K_Nsp.cdf</filename><filesize>194455495</filesize><checksum>59ae263311a2cf63b8d1b9b4cc7d663b</checksum></srcFile1><srcFile2><filename>Mapping250K_Nsp.na26.annot.csv</filename><filesize>447194757</filesize><checksum>5f651b351a0d97e5e9e657f123e6ed05</checksum></srcFile2></srcFiles>
##  Chip type: Mapping250K_Nsp
##  Platform: Affymetrix

##                  snp cnp affxSnp other  total
##  enzyme1-only 261563   0       0     0 261563
##  missing         701   0       0    74    775
##  total        262264   0       0    74 262338
