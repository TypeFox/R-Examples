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
# Import UGP from CSV files
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
tags <- sprintf("na%s,%s,%s%s", naVersion, genomeVersion, user, datestamp);
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
ugp <- AromaUgpFile$byChipType("GenomeWideSNP_6,Full", tags="na31");
ugp0 <- AromaUgpFile$byChipType("GenomeWideSNP_6,Full", tags="na30");

print(table(ugp[,1], exclude=NULL));
print(table(ugp0[,1], exclude=NULL));

## GenomeWideSNP_6,Full,na31,hg19,HB20110328
##      1      2      3      4      5      6      7      8
## 146626 153717 127767 120377 115735 112727 100925  98308
##      9     10     11     12     13     14     15     16
##  82300  93635  89598  87341  65927  57115  53596  54111
##     17     18     19     20     21     22     23     24
##  46609  52102  30365  43649  25105  24438  87271   9688
##     25   <NA>
##    455   1928

## GenomeWideSNP_6,na30,hg18,HB20100215.ugp
##      1      2      3      4      5      6      7      8
## 146401 153663 127766 120296 115672 112825 100996  98277
##      9     10     11     12     13     14     15     16
##  82168  93592  89525  87321  66067  57103  53556  54182
##     17     18     19     20     21     22     23     24
##  46632  52093  30299  43628  25111  24484  87200   9483
##     25   <NA>
##    445   2630


rr <- whichVector(ugp[,1] != ugp0[,1]);
str(rr);
## integer(0)

rr <- whichVector(ugp[,2] != ugp0[,2]);
str(rr);
## integer(0)
