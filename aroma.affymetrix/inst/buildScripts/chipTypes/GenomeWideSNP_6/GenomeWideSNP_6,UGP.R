if (interactive()) savehistory();
library("aroma.affymetrix");
library("R.menu");
verbose <- Verbose(threshold=-10, timestamp=TRUE);
options(width=60);

chipType <- "GenomeWideSNP_6";
cdfTags <- "Full";


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
naVersions <- 30:33;
naVersion <- textMenu(naVersions, title="Choose NetAffx version: ", value=TRUE);
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

genomeVersion <- getHeaderAttributes(csvList[[1]])[["genome-version"]];
genomeVersion <- Arguments$getCharacter(genomeVersion, nchar=4, length=c(1,1));
dbSNPVersion <- getHeaderAttributes(csvList[[1]])[["dbSNP-version"]];
dbSNPVersion <- Arguments$getInteger(dbSNPVersion);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Import UGP from CSV files
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
tags <- sprintf("na%s,%s,dbSNP%d,%s%s", naVersion, genomeVersion, dbSNPVersion, user, datestamp);
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
  units <- importFrom(ugp, csv, verbose=verbose);
  str(units);
  ## GenomeWideSNP_6.na<XX>.annot.csv:    int [1:934968] 334945 334944 ...
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

footer <- readFooter(ugp);
footer$createdBy <- list(
  fullname = fullname,
  email = email
);
names(srcFileTags) <- sprintf("srcFile%d", seq(along=srcFileTags));
footer$srcFiles <- srcFileTags;
writeFooter(ugp, footer);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Statistics
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
print(ugp);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# WHAT'S NEW:
#
# o na32 -> na33 (hg19, dbSNP 132 -> dbSNP 137)
#   No differences.
# o na31 -> na32 (hg19, dbSNP 131 -> dbSNP 132)
#   89 SNPs have jumped chromosomes. 422 SNPs have moved within
#   chromosomes.  This is most likely due to updates in dbSNP.
#   In addition there are several migrations from/to being NA.
# o na30 -> na31
#   Lots of changes, especially large changes in positions.
# o na27.1 -> na30
#   No differences
# o na27 -> na27.1
#   No differences
# o na26 -> na27
#   Two units (932039, 1872834) where moved from ChrX to ChrY.
#   Same location.
# o na24 -> na26
#   Only minor modifications for non-missing values:
#   - three loci changed chromosomes
#   - an additional 23 loci changed positions, of which only 17 moved
#     more than 2 base pairs.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ugp <- AromaUgpFile$byChipType("GenomeWideSNP_6,Full", tags="na33");
ugp0 <- AromaUgpFile$byChipType("GenomeWideSNP_6,Full", tags="na32");

print(table(ugp[,1], exclude=NULL));
print(table(ugp0[,1], exclude=NULL));

## GenomeWideSNP_6,Full,na33,hg19,dbSNP137,HB20140117
##      1      2      3      4      5      6      7      8
## 146458 153641 127776 120351 115709 112858 100828  98274
##      9     10     11     12     13     14     15     16
##  82178  93627  89604  87336  66091  57118  53554  54181
##     17     18     19     20     21     22     23     24
##  46590  52103  30362  43647  25091  24405  87102   9442
##     25   <NA>
##    411   2678

## GenomeWideSNP_6,Full,na32,hg19,dbSNP132,HB20140118
##      1      2      3      4      5      6      7      8
## 146458 153641 127776 120351 115709 112858 100828  98274
##      9     10     11     12     13     14     15     16
##  82178  93627  89604  87336  66091  57118  53554  54181
##     17     18     19     20     21     22     23     24
##  46590  52103  30362  43647  25091  24405  87102   9442
##     25   <NA>
##    411   2678

## GenomeWideSNP_6,Full,na31,hg19,HB20110328
##      1      2      3      4      5      6      7      8
## 146626 153717 127767 120377 115735 112727 100925  98308
##      9     10     11     12     13     14     15     16
##  82300  93635  89598  87341  65927  57115  53596  54111
##     17     18     19     20     21     22     23     24
##  46609  52102  30365  43649  25105  24438  87271   9688
##     25   <NA>
##    455   1928



units <- whichVector(ugp[,1] != ugp0[,1]);
str(units);
## int(0)

units <- whichVector(ugp[,2] != ugp0[,2]);
str(units);
## int(0)
