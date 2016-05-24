source("init.R");

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Settings
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
naVersion <- "25";
user <- "HB";
datestamp <- "20080424";

chipType <- "GenomeWideSNP_5";
cdfTags <- "Full,r2";
nbrOfEnzymes <- 2;


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup required annotation files
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if (!exists("cdf")) {
  cdf <- AffymetrixCdfFile$byChipType(chipType, tags=cdfTags);
  rm(csvList);
}
print(cdf);

if (!exists("csvList", mode="list")) {
  csvList <- list();
  tags <- sprintf(".na%s", naVersion);
  csv <- AffymetrixNetAffxCsvFile$byChipType(chipType, tags=tags);
  csvList[["main"]] <- csv;
  tags <- sprintf(".cn.na%s", naVersion);
  csv <- AffymetrixNetAffxCsvFile$byChipType(chipType, tags=tags);
  csvList[["CN"]] <- csv;
  rm(csv, tags);
}
print(csvList);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Import UFL from CSV files
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
tags <- sprintf("na%s,%s%s", naVersion, user, datestamp);
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
  units <- importFrom(ufl, csv, verbose=log);
  str(units);
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Statistics
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
x <- summaryOfUnits(ufl, enzymes=c("NspI", "StyI"));
print(x);
##                  snp    cnp affxSnp other  total
##  enzyme1-only 116979 140099       0     0 257078
##  enzyme2-only  74135   1208       0     0  75343
##  both         248980 171077       0     0 420057
##  missing       60474 104885    3022    69 168450
##  total        500568 417269    3022    69 920928

cdf2 <- AffymetrixCdfFile$byChipType("GenomeWideSNP_6", tags="Full");
tags <- sprintf("na%s,%s%s", naVersion, user);
ufl2 <- AromaUflFile$byChipType(getChipType(cdf2), tags=tags);

units2 <- match(getUnitNames(cdf), getUnitNames(cdf2));
keep <- is.finite(units2);
print(summary(keep));
units <- which(keep);
units2 <- units2[keep];
for (ee in 1:nbrOfEnzymes(ufl))
 ufl[units,ee] <- ufl2[units2,ee];

##                 snp    cnp affxSnp other  total
## enzyme1-only 129848 148340       0     0 278188
## enzyme2-only  90353   1208       0     0  91561
## both         273682 181097       0     0 454779
## missing        6685  86624    3022    69  96400
## total        500568 417269    3022    69 920928



