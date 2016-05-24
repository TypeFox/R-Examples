library("aroma.affymetrix");
verbose <- Verbose(threshold=-10, timestamp=TRUE);
options(width=60);

chipType <- "GenomeWideSNP_6";
cdfTags <- "Full";
nbrOfEnzymes <- 2;


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# User settings
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
fullname <- "Henrik Bengtsson";
user <- "HB";
email <- getOption(aromaSettings, "user/email");
if (is.null(email)) {
  email <- "JohnDoe@email.com";
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Settings
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
naVersion <- "32";
genomeVersion <- "hg19";
datestamp <- format(Sys.Date(), format="%Y%m%d");



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Setup required annotation files
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
cdf <- AffymetrixCdfFile$byChipType(chipType, tags=cdfTags);
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
tags <- c("TEST-ONLY", tags);
ufl <- AromaUflFile$allocateFromCdf(cdf, tags=tags, nbrOfEnzymes=nbrOfEnzymes, overwrite=TRUE);
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


# Try to load the created annotation file
ufl <- AromaUflFile$byChipType("GenomeWideSNP_6,Full", tags="TEST-ONLY");
print(ufl);
x <- summaryOfUnits(ufl);
print(x);
