if (interactive()) savehistory();
library("aroma.affymetrix");
library("R.menu");
verbose <- Verbose(threshold=-10, timestamp=TRUE);
options(width=60);

chipType <- "MoGene-1_0-st-v1";


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
genomeVersions <- c("32"="mm8");
naVersions <- names(genomeVersions);
choices <- sprintf("na%s (%s)", naVersions, genomeVersions);
choice <- textMenu(choices, title="Choose NetAffx version: ", value=FALSE);
naVersion <- naVersions[choice];
genomeVersion <- genomeVersions[naVersion];
datestamp <- format(Sys.Date(), format="%Y%m%d");



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Setup required annotation files
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
cdf <- AffymetrixCdfFile$byChipType(chipType);
print(cdf);

tags <- sprintf(".na%s", naVersion);

pattern <- sprintf("^%s%s.*[.]transcript.csv$", chipType, tags);
pathname <- AffymetrixNetAffxCsvFile$findByChipType(chipType, tags=tags, pattern=pattern);
if (isFile(pathname)) {
  csv <- AffymetrixNetAffxCsvFile(pathname);
}
rm(tags);
print(csv);


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

# Read unit names

colClasses <- c("*"="NULL", "probesetId"="character", "seqname"="character", "(start|stop)"="character");
verbose && enter(verbose, "Read annotation data for all units");
data <- readDataFrame(csv, colClasses=colClasses);
# Coerce to integers
data$start <- as.integer(data$start);
data$stop <- as.integer(data$stop);
verbose && str(verbose, data);
verbose && exit(verbose);

verbose && enter(verbose, "Map to unit names");
units <- indexOf(cdf, names=data$probesetId);
verbose && print(verbose, summary(units));
verbose && str(verbose, units);
verbose && exit(verbose);

verbose && enter(verbose, "Drop unknown units");
keep <- is.finite(units);
units <- units[keep];
data <- data[keep,];
verbose && cat(verbose, "Number of known units: ", length(units));
verbose && str(verbose, data);
verbose && exit(verbose);

verbose && enter(verbose, "Parse annotation data");
chr <- data$seqname;
chr <- gsub("chr", "", chr, fixed=TRUE);
chr <- gsub("_random", "", chr, fixed=TRUE);
chr[chr == "X"] <- 23L;
chr[chr == "Y"] <- 24L;
chr[chr == "M"] <- 25L;
chr <- as.integer(chr);
position <- (data$start + data$stop) / 2;
verbose && print(verbose, table(chr));
verbose && exit(verbose);


verbose && enter(verbose, "Write to file");
ugp[units,1] <- chr;
ugp[units,2] <- position;
verbose && exit(verbose);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Update the file footer
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
srcFileTags <- list();
srcFiles <- c(list(cdf), list(csv));
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
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
ugp <- AromaUgpFile$byChipType("MoGene-1_0-st-v1,r3", tags="na32");
ugp0 <- AromaUgpFile$byChipType("MoGene-1_0-st-v1,r3", tags="na32");

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
