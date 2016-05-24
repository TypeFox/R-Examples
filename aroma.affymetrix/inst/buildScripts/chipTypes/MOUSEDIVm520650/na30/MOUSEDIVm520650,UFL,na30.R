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
nbrOfEnzymes <- 2;

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
# Setup/Allocate UFL
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
tags <- sprintf("na%s,%s,%s%s", naVersion, genomeVersion, user, datestamp);
ufl <- NULL;
tryCatch({
  ufl <- AromaUflFile$byChipType(getChipType(cdf), tags=tags);
}, error = function(ex) {})
if (is.null(ufl)) {
  ufl <- AromaUflFile$allocateFromCdf(cdf, tags=tags, nbrOfEnzymes=nbrOfEnzymes);
}
print(ufl);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Import from CSV files
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
for (kk in seq(along=csvList)) {
  csv <- csvList[[kk]];
  print(csv);
  units <- importFrom(ufl, csv, verbose=verbose);
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


footer <- readFooter(ufl);
footer$createdBy <- list(
  fullname = "Henrik Bengtsson", 
  email = sprintf("%s@%s", "henrik.bengtsson", "aroma-project.org")
);
names(srcFileTags) <- sprintf("srcFile%d", seq(along=srcFileTags));
footer$srcFiles <- srcFileTags;
writeFooter(ufl, footer);

print(ufl);



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# WHAT'S NEW:
#
# o ... -> na30
#   - First one for this chip type.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
ufl <- AromaUflFile$byChipType("MOUSEDIVm520650", tags="na30");
x <- summaryOfUnits(ufl);

print(x);
##              snp cnp affxSnp   other   total
## enzyme1-only   0   0       0       0       0
## enzyme2-only   0   0       0       0       0
## both           0   0       0 2227922 2227922
## missing        0   0       0  230775  230775
## total          0   0       0 2458697 2458697

## TODO: summaryOfUnits() should utilize getUnitTypes(cdf).
