library("aroma.affymetrix");
log <- Arguments$getVerbose(-4, timestamp=TRUE);



dataSet <- "HapMap,CEU,testset";
chipTypes <- c("Mapping50K_Hind240", "Mapping50K_Xba240");

# Expected sample names
sampleNames <- c("NA06985", "NA06991", "NA06993", 
                 "NA06994", "NA07000", "NA07019");

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Tests for setting up CEL sets and locating the CDF file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
csRList <- list();
for (chipType in chipTypes) {
  cs <- AffymetrixCelSet$byName(dataSet, chipType=chipType);
  print(cs);
  stopifnot(identical(getNames(cs), sampleNames));
  csRList[[chipType]] <- cs;
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Probe-level modelling test #1
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Do it on the raw data so that we can validate the results
csList <- csRList;  

plmList <- list();
for (chipType in names(csList)) {
  cs <- csList[[chipType]];
  plm <- RmaPlm(cs);
  print(plm);
  fit(plm, verbose=log);
  ces <- getChipEffectSet(plm);
  print(ces);
  stopifnot(identical(getNames(ces), getNames(cs)));

  # Assert a few results
  if (chipType == "Mapping250K_Nsp") {
    unit <- 100;
    theta <- readUnits(ces, units=unit, verbose=log);
    data <- unlist(theta, use.names=FALSE);
    data <- sum(data, na.rm=TRUE);
    if (!isZero(data - 57448.85473633, eps=.Machine$float.eps)) {
      print(theta);
      throw("Non-reproducible chip-effect estimates detected for unit ", unit);
    }
  }

  plmList[[chipType]] <- plm;
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# PLM residual test
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
rsList <- list();
for (chipType in names(plmList)) {
  plm <- plmList[[chipType]];
  print(plm);
  rs <- calculateResidualSet(plm, verbose=log);
  print(rs);
  stopifnot(identical(getNames(rs), getNames(cs)));

  # Assert a few results
  if (chipType == "Mapping50K_Hind240") {
    cs <- getDataSet(plm);
    ces <- getChipEffectSet(plm);
    paf <- getProbeAffinityFile(plm);
    unit <- 500;
    y <- readUnits(plm, units=unit, verbose=log);
    theta <- readUnits(ces, units=unit, verbose=log);
    phi <- readUnits(paf, units=unit, verbose=log);
    resids <- readUnits(rs, units=unit, verbose=log);
    str(resids);
  }

  rsList[[chipType]] <- rs;
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Spatial residual plots test
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ae <- ArrayExplorer(rsList);
setColorMaps(ae, c("log2,log2neg,rainbow", "log2,log2pos,rainbow"));
print(ae);
stopifnot(identical(unname(getArrays(ae)), getNames(cs)));
process(ae, arrays=1:2, interleaved="auto", verbose=log);



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# PLM weights (from residuals) test
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
for (chipType in names(plmList)[1]) {
  plm <- plmList[[chipType]];
  print(plm);
  rs <- getResidualSet(plm, verbose=log);
  print(rs);
  stopifnot(identical(getNames(rs), getNames(cs)));

  ws <- calculateWeights(plm, verbose=log);
  print(rs);

  rsList[[chipType]] <- rs;
}
