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

ptb <- AffymetrixProbeTabFile$byChipType(chipType);
print(ptb);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Import UGP from PTB files
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
tags <- sprintf("na%s,%s,%s%s", naVersion, genomeVersion, user, datestamp);
tags <- sprintf("%s%s", user, datestamp);
tags <- c("build36", "mm8", tags);
ugp <- NULL;
tryCatch({
  acp <- AromaCellPositionFile$byChipType(getChipType(cdf), tags=tags);
}, error = function(ex) {})
if (is.null(ugp)) {
  acp <- AromaCellPositionFile$allocateFromCdf(cdf, tags=tags);
}
print(acp);

# Read unit names

colClasses <- c("*"="NULL", "probeID"="integer", "seqname"="character", "(start|stop)"="integer");
verbose && enter(verbose, "Read annotation data for all units");
data <- readDataFrame(ptb, colClasses=colClasses);
verbose && str(verbose, data);
verbose && exit(verbose);

verbose && enter(verbose, "Parse annotation data");
cells <- data$probeID;
cells <- Arguments$getIndices(cells, max=nbrOfCells(cdf));
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
acp[cells,1] <- chr;
acp[cells,2] <- position;
verbose && exit(verbose);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Update the file footer
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
srcFileTags <- list();
srcFiles <- c(list(cdf), list(ptb));
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

footer <- readFooter(acp);
footer$createdBy <- list(
  fullname = fullname, 
  email = email
);
names(srcFileTags) <- sprintf("srcFile%d", seq(along=srcFileTags));
footer$srcFiles <- srcFileTags;
writeFooter(acp, footer);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Statistics
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
print(acp);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# WHAT'S NEW:
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
acp <- AromaUgpFile$byChipType("MoGene-1_0-st-v1,r3");
acp0 <- AromaAcpFile$byChipType("MoGene-1_0-st-v1,r3");

print(table(acp[,1], exclude=NULL));
print(table(acp0[,1], exclude=NULL));

## GenomeWideSNP_6,Full,na31,hg19,HB20110328
##      1      2      3      4      5      6      7      8
## 146626 153717 127767 120377 115735 112727 100925  98308
##      9     10     11     12     13     14     15     16
##  82300  93635  89598  87341  65927  57115  53596  54111
##     17     18     19     20     21     22     23     24
##  46609  52102  30365  43649  25105  24438  87271   9688
##     25   <NA>
##    455   1928

## GenomeWideSNP_6,na30,hg18,HB20100215.acp
##      1      2      3      4      5      6      7      8
## 146401 153663 127766 120296 115672 112825 100996  98277
##      9     10     11     12     13     14     15     16
##  82168  93592  89525  87321  66067  57103  53556  54182
##     17     18     19     20     21     22     23     24
##  46632  52093  30299  43628  25111  24484  87200   9483
##     25   <NA>
##    445   2630


rr <- whichVector(acp[,1] != acp0[,1]);
str(rr);
## integer(0)

rr <- whichVector(acp[,2] != acp0[,2]);
str(rr);
## integer(0)
