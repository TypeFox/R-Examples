if (interactive()) savehistory();
library("aroma.affymetrix");
library("R.menu");
verbose <- Verbose(threshold=-10, timestamp=TRUE);
options(width=60);

chipType <- "GenomeWideSNP_6";


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# User settings
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## setOption(aromaSettings, "user/initials", "HB");
## setOption(aromaSettings, "user/fullname", "Henrik Bengtsson");
## obf <- sprintf("%s@%s", "henrik.bengtsson", "aroma-project.org");
## setOption(aromaSettings, "user/email", obf);
## saveAnywhere(aromaSettings);

fullname <- getOption(aromaSettings, "user/fullname");
stopifnot(!is.null(fullname));
email <- getOption(aromaSettings, "user/email");
stopifnot(!is.null(email));
user <- getOption(aromaSettings, "user/initials");
stopifnot(!is.null(user));


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup required annotation files
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cdfD <- AffymetrixCdfFile$byChipType(chipType);
print(cdfD);

cdf <- AffymetrixCdfFile$byChipType(chipType, tags="Full");
print(cdf);

unitNamesD <- getUnitNames(cdfD);
unitMap <- indexOf(cdf, names=unitNamesD);
# Sanity check
stopifnot(!anyMissing(unitMap));

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Allocate UFL and UGP for the default CDF
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
naTag <- sprintf("na%s", naVersion);
tags <- sprintf("%s,%s,dbSNP%d,%s%s", naTag, genomeVersion, dbSNPVersion, user, datestamp);

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
  footer$createdBy = list(
    fullname = fullname,
    email = email
  );
  keep <- (names(srcFileTags) %in% c("cdfD", "cdf", what));
  srcFileTagsT <- srcFileTags[keep];
  names(srcFileTagsT) <- sprintf("srcFile%d", seq(along=srcFileTagsT));
  footer$srcFiles <- srcFileTagsT;
  writeFooter(afD, footer);
}
