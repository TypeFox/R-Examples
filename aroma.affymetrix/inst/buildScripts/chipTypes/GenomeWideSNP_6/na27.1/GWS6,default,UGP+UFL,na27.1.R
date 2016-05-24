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

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Setup required annotation files
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
if (!exists("cdfD")) {
  cdfD <- AffymetrixCdfFile$byChipType(chipType);
}
print(cdfD);

if (!exists("cdf")) {
  cdf <- AffymetrixCdfFile$byChipType(chipType, tags="Full");
}
print(cdf);


unitNamesD <- getUnitNames(cdfD);
unitMap <- indexOf(cdf, names=unitNamesD);
# Sanity check
stopifnot(!anyMissing(unitMap));

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Allocate UFL and UGP for the default CDF
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
naTag <- sprintf("na%s", naVersion);
tags <- sprintf("%s,%s%s", naTag, user, datestamp);

ufl <- AromaUflFile$byChipType(getChipType(cdf), tags=naTag);
uflD <- AromaUflFile$allocateFromCdf(cdfD, tags=tags, nbrOfEnzymes=nbrOfEnzymes(ufl));
print(uflD);

ugp <- AromaUgpFile$byChipType(getChipType(cdf), tags=naTag);
ugpD <- AromaUgpFile$allocateFromCdf(cdfD, tags=tags);
print(ugpD);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Import data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
for (what in c("ufl", "ugp")) {
  af <- get(what);
  afD <- get(sprintf("%sD", what));
  for (cc in seq(length=nbrOfColumns(afD))) {
    afD[,cc] <- af[unitMap,cc];
  }
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Update the file footer
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
if (!exists("srcFileTags", mode="list")) {
  srcFiles <- list(cdfD=cdfD, cdf=cdf, ufl=ufl, ugp=ugp);
  srcFileTags <- lapply(srcFiles, FUN=function(srcFile) {
    list(
      filename=getFilename(srcFile), 
      filesize=getFileSize(srcFile), 
      checksum=getChecksum(srcFile)
    );
  });
  print(srcFileTags);
}

for (what in c("ufl", "ugp")) {
  af <- get(what);
  afD <- get(sprintf("%sD", what));
  footer <- readFooter(afD);
  footer$createdOn <- format(Sys.time(), "%Y%m%d %H:%M:%S", usetz=TRUE);
  footer$createdBy = list(
    fullname = "Henrik Bengtsson", 
    email = sprintf("%s@%s", "henrik.bengtsson", "aroma-project.org")
  );
  keep <- (names(srcFileTags) %in% c("cdfD", "cdf", what));
  srcFileTagsT <- srcFileTags[keep];
  names(srcFileTagsT) <- sprintf("srcFile%d", seq(along=srcFileTagsT));
  footer$srcFiles <- srcFileTagsT;
  writeFooter(afD, footer);
}
