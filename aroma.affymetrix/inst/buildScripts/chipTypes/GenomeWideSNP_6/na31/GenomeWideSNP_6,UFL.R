if (interactive()) savehistory();
library("aroma.affymetrix");
library("R.menu");
verbose <- Verbose(threshold=-10, timestamp=TRUE);
options(width=60);

chipType <- "GenomeWideSNP_6";
cdfTags <- "Full";
nbrOfEnzymes <- 2;


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
# Settings
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
genomeVersions <- c("30"="hg18", "31"="hg19");
naVersions <- names(genomeVersions);
choices <- sprintf("na%s (%s)", naVersions, genomeVersions);
choice <- textMenu(choices, title="Choose NetAffx version: ", value=FALSE);
naVersion <- naVersions[choice];
genomeVersion <- genomeVersions[naVersion];
datestamp <- format(Sys.Date(), format="%Y%m%d");



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Setup required annotation files
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
cdf <- AffymetrixCdfFile$byChipType(chipType, tags=cdfTags);
rm(csvList);
print(cdf);

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
print(csvList);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Import UFL from CSV files
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

for (kk in seq(along=csvList)) {
  csv <- csvList[[kk]];
  print(csv);
  units <- importFrom(ufl, csv, verbose=verbose);
  str(units);
  ## GenomeWideSNP_6.na<XX>.annot.csv:  int [1:934968] 334945 334944 ...
  ## GenomeWideSNP_6.cn.na<XX>.annot.csv: int [1:945826] 935622 935777 ...
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Update the file footer
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
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

footer <- readFooter(ufl);
footer$createdBy <- list(
  fullname = fullname, 
  email = email
);
names(srcFileTags) <- sprintf("srcFile%d", seq(along=srcFileTags));
footer$srcFiles <- srcFileTags;
writeFooter(ufl, footer);

print(ufl);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# WHAT'S NEW:
#
# o na30 -> na31
#   - Lots of differences; Affymetrix no longer filter for us.
# o na27.1 -> na30
#   - No changes.
# o na27 -> na27.1
#   - No changes.
# o na26 -> na27
#   - No changes.
# o na24 -> na26
#   - All changes are for SNP units.
#   - There are 6 NspI and 14 StyI changes in SNP fragment lengths,
#     which some are only minor.
#   - There are 1108 more SNPs that now have missing values.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
ufl <- AromaUflFile$byChipType("GenomeWideSNP_6,Full", tags="na31");
ufl0 <- AromaUflFile$byChipType("GenomeWideSNP_6,Full", tags="na30");
x <- summaryOfUnits(ufl);
x0 <- summaryOfUnits(ufl0);

print(x);
print(x0);

## GenomeWideSNP_6,na31,hg19,HB20110328.ufl
##                 snp    cnp affxSnp other   total
## enzyme1-only    577   9609       0     0   10186
## enzyme2-only    502  23432       0     0   23934
## both         929322 912574       0     0 1841896
## missing        1545    211    3022   621    5399
## total        931946 945826    3022   621 1881415

## GenomeWideSNP_6,na30,hg18,HB20100215.ufl
##                  snp    cnp affxSnp other   total
##  enzyme1-only 246080 451191       0     0  697271
##  enzyme2-only 160899      0       0     0  160899
##  both         522472 494615       0     0 1017087
##  missing        2495     20    3022   621    6158
##  total        931946 945826    3022   621 1881415

## GenomeWideSNP_6,na27.1,hg??,HB20??????.ufl
##                  snp    cnp affxSnp other   total
##  enzyme1-only 246080 451191       0     0  697271
##  enzyme2-only 160899      0       0     0  160899
##  both         522472 494615       0     0 1017087
##  missing        2495     20    3022   621    6158
##  total        931946 945826    3022   621 1881415

print(x-x0);
##               snp cnp affxSnp other total
##  enzyme1-only   0   0       0     0     0
##  enzyme2-only   0   0       0     0     0
##  both           0   0       0     0     0
##  missing        0   0       0     0     0
##  total          0   0       0     0     0

# Differences
for (cc in 1:nbrOfColumns(ufl)) {
  units <- whichVector(ufl[,cc] != ufl0[,cc], na.rm=TRUE);
  if (length(units) > 0) {
    df <- cbind(units, ufl[units,cc], ufl0[units,cc], 
                       ufl[units,cc]-ufl0[units,cc]);
    colnames(df)[ncol(df)] <- "delta";
    print(df);
  }
}

# Differences in NAs
for (cc in 1:nbrOfColumns(ufl)) {
  units <- whichVector(is.na(ufl[,cc]) != is.na(ufl0[,cc]));
  if (length(units) > 0) {
    df <- cbind(units, ufl[units,cc], ufl0[units,cc]);
    str(df);
  }
}
