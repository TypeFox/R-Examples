library("aroma.core");
log <- Arguments$getVerbose(-20, timestamp=TRUE);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Allocate ACF file set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dataSet <- "FooDataSet";
platform <- "Generic";
chipType <- "FooBar";
nbrOfRows <- 402;
path <- filePath("callData", dataSet, chipType);

for (kk in 1:6) {
  filename <- sprintf("foo%02d,genotypes.acf", kk);
  acf <- AromaUnitGenotypeCallFile$allocate(filename, path=path, nbrOfRows=nbrOfRows, platform=platform, chipType=chipType, overwrite=TRUE);
  acf2 <- AromaUnitGenotypeCallFile(filename, path=path);
  stopifnot(equals(acf2, acf));
}

acs <- AromaUnitGenotypeCallSet$byName(dataSet, chipType=chipType);
print(acs);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Updating and extract genotypes
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Calls by user
calls <- c(NA, "NA", "NC", "", "-", "A", "B", "AA", "AB", "BB", "AAA", "AAB", "BAAB");

# Truth to compare with
callsT <- calls;
callsT <- strsplit(callsT, split="", fixed=TRUE);
callsT <- lapply(callsT, FUN=sort);
callsT <- sapply(callsT, FUN=paste, collapse="");
callsT[calls == "-"] <- "";
callsT[calls == "NC"] <- "NC";
callsT[(is.na(calls) | calls == "NA")] <- NA;

units <- seq(from=1, length=length(calls), by=2);
for (kk in seq_along(acs)) {
  acf <- acs[[kk]];
  updateGenotypes(acf, units=units, calls=calls);
  calls2 <- extractGenotypes(acf, units=units, drop=TRUE);
  stopifnot(identical(calls2, callsT));
  calls <- c(calls[-1], calls[1]);
  callsT <- c(callsT[-1], callsT[1]);
}

calls <- extractGenotypes(acs, units=units);
print(calls);
