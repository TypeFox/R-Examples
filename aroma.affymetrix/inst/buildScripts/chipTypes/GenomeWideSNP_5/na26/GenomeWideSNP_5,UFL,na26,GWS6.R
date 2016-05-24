if (interactive()) savehistory();
library("aroma.affymetrix");
log <- Verbose(threshold=-10, timestamp=TRUE);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Settings
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
naVersion <- "26";
user <- "HB";
datestamp <- "20080723";

chipType <- "GenomeWideSNP_5";
cdfTags <- "Full,r2";
nbrOfEnzymes <- 2;

chipType6 <- "GenomeWideSNP_6,Full";

footer <- list(
  createdOn = format(Sys.time(), "%Y%m%d %H:%M:%S", usetz=TRUE),
  createdBy = list(
    fullname = "Henrik Bengtsson", 
    email = sprintf("%s@%s", "henrik.bengtsson", "aroma-project.org")
  ),
  srcFiles = list()
);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Setup required annotation files
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
if (!exists("cdf")) {
  cdf <- AffymetrixCdfFile$byChipType(chipType, tags=cdfTags);
  rm(csvList);
}
print(cdf);

if (!exists("cdf6")) {
  cdf6 <- AffymetrixCdfFile$byChipType(chipType6);
}
print(cdf6);

if (!exists("ufl6")) {
  ufl6Tags <- sprintf("na%s", naVersion);
  ufl6 <- AromaUflFile$byChipType(getChipType(cdf6), tags=ufl6Tags);
}
print(ufl6);



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Setup UFL file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
tags <- sprintf("na%s+,%s%s", naVersion, user, datestamp);
if (!exists("ufl")) {
  ufl <- AromaUflFile$byChipType(getChipType(cdf), tags=tags);
}
print(ufl);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Import from the GWS6 file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Identify units with missing UFL data
isMissing <- (is.na(ufl[,1]) & is.na(ufl[,2]));
print(summary(drop(isMissing)));
##    Mode   FALSE    TRUE
## logical  752478  168450

units <- which(isMissing);
unitNames <- getUnitNames(cdf, units=units);

# Identify which exists in GWS6
units6 <- match(unitNames, getUnitNames(cdf6));
keep <- is.finite(units6);
print(summary(keep));
##    Mode   FALSE    TRUE
## logical   92097   76353
units <- which(keep);
units6 <- units6[keep];
rm(keep);

# Import data
for (ee in 1:nbrOfEnzymes(ufl)) {
 ufl[units,ee] <- ufl6[units6,ee];
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Update footer
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
footer <- readFooter(ufl);
srcFile <- list(
  filename=getFilename(ufl6), 
  filesize=getFileSize(ufl6), 
  checksum=getChecksum(ufl6)
);
footer$srcFiles <- c(footer$srcFiles, list(srcFile=srcFile));
writeFooter(ufl, footer);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Statistics
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
print(ufl);
## AromaUflFile:
## Name: GenomeWideSNP_5
## Tags: Full,r2,na26,HB20080723
## Pathname: annotationData/chipTypes/GenomeWideSNP_5/GenomeWideSNP_5,Full,r2,na26,HB20080723.ufl
## File size: 3.51MB
## RAM: 0.00MB
## Number of data rows: 920928
## File format: v1
## Dimensions: 920928x2
## Column classes: integer, integer
## Number of bytes per column: 2, 2
## Footer: <platform>Affymetrix</platform><chipType>GenomeWideSNP_5,Full,r2</chipType><createdOn>20080723 16:09:57 PDT</createdOn><createdBy><fullname>Henrik Bengtsson</fullname><email>[...]</email></createdBy><srcFiles><srcFile><filename>GenomeWideSNP_5.na26.annot.csv</filename><filesize>755337946</filesize><checksum>af59235b6fccada7f871257149a89215</checksum></srcFile><srcFile><filename>GenomeWideSNP_5.na26.annot.csv</filename><filesize>755337946</filesize><checksum>af59235b6fccada7f871257149a89215</checksum></srcFile><srcFile><filename>GenomeWideSNP_5.na26.annot.csv</filename><filesize>755337946</filesize><checksum>af59235b6fccada7f871257149a89215</checksum></srcFile></srcFiles>
## Chip type: GenomeWideSNP_5,Full,r2
## Platform: Affymetrix
print(getChecksum(ufl));
## [1] "9d73f5a7689cac14c331ab661cf50ea9"

x <- summaryOfUnits(ufl, enzymes=c("NspI", "StyI"));
print(x);
## GenomeWideSNP_5,Full,na26+,HB20080723:
##                 snp    cnp affxSnp other  total
## enzyme1-only 102627 140099       0     0 242726
## enzyme2-only  90353   1208       0     0  91561
## both         249553 171077       0     0 420630
## missing       58035 104885    3022    69 166011
## total        500568 417269    3022    69 920928
