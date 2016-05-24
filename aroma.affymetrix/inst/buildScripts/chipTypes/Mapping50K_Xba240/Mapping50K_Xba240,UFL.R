if (interactive()) savehistory();
library("aroma.affymetrix");
verbose <- Verbose(threshold=-50, timestamp=TRUE);
options(width=60);


chipType <- "Mapping50K_Xba240";


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# User settings
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
## setOption(aromaSettings, "user/initials", "HB");
## setOption(aromaSettings, "user/fullname", "Henrik Bengtsson");
## email <- sprintf("%s@%s", "henrik.bengtsson", "aroma-project.org")
## setOption(aromaSettings, "user/email", email);
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


enzyme <- c("Mapping50K_Hind240"="HindIII", "Mapping50K_Xba240"="XbaI")[chipType];
print(enzyme);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Setup required annotation files
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
cdf <- AffymetrixCdfFile$byChipType(chipType);
rm(csv);
print(cdf);

tags <- sprintf(".na%s", naVersion);
pathname <- AffymetrixNetAffxCsvFile$findByChipType(chipType, tags=tags);
if (isFile(pathname)) {
  csv <- AffymetrixNetAffxCsvFile(pathname);
}
rm(tags);
print(csv);



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Import UFL from CSV files
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
tags <- sprintf("na%s,%s,%s%s", naVersion, genomeVersion, user, datestamp);
ufl <- NULL;
tryCatch({
  ufl <- AromaUflFile$byChipType(getChipType(cdf), tags=tags);
}, error = function(ex) {})
if (is.null(ufl)) {
  ufl <- AromaUflFile$allocateFromCdf(cdf, tags=tags);
}
print(ufl);


stopifnot(!is.na(enzyme));
units <- importFrom(ufl, csv, enzymes=enzyme, verbose=verbose);
str(units);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Update the file footer
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
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

footer <- readFooter(ufl);
footer$createdOn <- format(Sys.time(), "%Y%m%d %H:%M:%S", usetz=TRUE);
footer$createdBy = list(
  fullname = fullname, 
  email = email
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

## AromaUflFile:
## Name: Mapping50K_Xba240
## Tags: na31,hg19,HB20110328
## Full name: Mapping50K_Xba240,na31,hg19,HB20110328
## Pathname: annotationData/chipTypes/Mapping50K_Xba240/Mapping50K_Xba240,na31,hg19
## ,HB20110328.ufl
## File size: 115.87 kB (118649 bytes)
## RAM: 0.01 MB
## Number of data rows: 59015
## File format: v1
## Dimensions: 59015x1
## Column classes: integer
## Number of bytes per column: 2
## Footer: <createdOn>20110328 19:03:23 PDT</createdOn><platform>Affymetrix</platfo
## rm><chipType>Mapping50K_Xba240</chipType><createdBy><fullname>Henrik Bengtsson</
## fullname><email>[...]</email></createdBy><srcFiles><srcFile1><fil
## ename>Mapping50K_Xba240.CDF</filename><filesize>57703961</filesize><checksum>7f0
## 26f6472f7721255717fb3f453be61</checksum></srcFile1><srcFile2><filename>Mapping50
## K_Xba240.na31.annot.csv</filename><filesize>84705402</filesize><checksum>4cc4997
## 35dd13177331994cde4e9f323</checksum></srcFile2></srcFiles>
## Chip type: Mapping50K_Xba240
## Platform: Affymetrix
## 
##                snp cnp affxSnp other total
## enzyme1-only 58530   0       0     0 58530
## missing        430   0       0    55   485
## total        58960   0       0    55 59015