# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Imports from affxparser
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
.invertMap <- function(...) {
  requireNamespace("affxparser") || throw("Package not loaded: affxparser")
  affxparser::invertMap(...)
}

.applyCdfGroups <- function(...) {
  requireNamespace("affxparser") || throw("Package not loaded: affxparser")
  affxparser::applyCdfGroups(...)
}

.cdfAddPlasqTypes <- function(...) {
  requireNamespace("affxparser") || throw("Package not loaded: affxparser")
  affxparser::cdfAddPlasqTypes(...)
}

.readBpmap <- function(...) {
  requireNamespace("affxparser") || throw("Package not loaded: affxparser")
  affxparser::readBpmap(...)
}

.readBpmapHeader <- function(...) {
  requireNamespace("affxparser") || throw("Package not loaded: affxparser")
  affxparser::readBpmapHeader(...)
}

.readCdfGroupNames <- function(...) {
  requireNamespace("affxparser") || throw("Package not loaded: affxparser")
  affxparser::readCdfGroupNames(...)
}

.readCdfHeader <- function(...) {
  requireNamespace("affxparser") || throw("Package not loaded: affxparser")
  affxparser::readCdfHeader(...)
}

.readCdf <- function(...) {
  requireNamespace("affxparser") || throw("Package not loaded: affxparser")
  affxparser::readCdf(...)
}

.readCdfQc <- function(...) {
  requireNamespace("affxparser") || throw("Package not loaded: affxparser")
  affxparser::readCdfQc(...)
}

.readCdfCellIndices <- function(...) {
  requireNamespace("affxparser") || throw("Package not loaded: affxparser")
  affxparser::readCdfCellIndices(...)
}

.readCdfUnits <- function(...) {
  requireNamespace("affxparser") || throw("Package not loaded: affxparser")
  affxparser::readCdfUnits(...)
}

.readCdfUnitNames <- function(...) {
  requireNamespace("affxparser") || throw("Package not loaded: affxparser")
  affxparser::readCdfUnitNames(...)
}

.readCdfIsPm <- function(...) {
  requireNamespace("affxparser") || throw("Package not loaded: affxparser")
  affxparser::readCdfIsPm(...)
}

.readCdfNbrOfCellsPerUnitGroup <- function(...) {
  requireNamespace("affxparser") || throw("Package not loaded: affxparser")
  affxparser::readCdfNbrOfCellsPerUnitGroup(...)
}

.readCelHeader <- function(...) {
  requireNamespace("affxparser") || throw("Package not loaded: affxparser")
  affxparser::readCelHeader(...)
}

.readCel <- function(...) {
  requireNamespace("affxparser") || throw("Package not loaded: affxparser")
  affxparser::readCel(...)
}

.readCelUnits <- function(...) {
  requireNamespace("affxparser") || throw("Package not loaded: affxparser")
  affxparser::readCelUnits(...)
}

.createCel <- function(...) {
  requireNamespace("affxparser") || throw("Package not loaded: affxparser")
  affxparser::createCel(...)
}

.compareCdfs <- function(...) {
  requireNamespace("affxparser") || throw("Package not loaded: affxparser")
  affxparser::compareCdfs(...)
}

.cdfHeaderToCelHeader <- function(...) {
  requireNamespace("affxparser") || throw("Package not loaded: affxparser")
  affxparser::cdfHeaderToCelHeader(...)
}

.findCdf <- function(...) {
  requireNamespace("affxparser") || throw("Package not loaded: affxparser")
  affxparser::findCdf(...)
}

.convertCdf <- function(...) {
  requireNamespace("affxparser") || throw("Package not loaded: affxparser")
  affxparser::convertCdf(...)
}

.convertCel <- function(...) {
  requireNamespace("affxparser") || throw("Package not loaded: affxparser")
  affxparser::convertCel(...)
}

.updateCel <- function(...) {
  requireNamespace("affxparser") || throw("Package not loaded: affxparser")
  affxparser::updateCel(...)
}

.updateCelUnits <- function(...) {
  requireNamespace("affxparser") || throw("Package not loaded: affxparser")
  affxparser::updateCelUnits(...)
}

.readCelRectangle <- function(...) {
  requireNamespace("affxparser") || throw("Package not loaded: affxparser")
  affxparser::readCelRectangle(...)
}

.readPgfEnv <- function(...) {
  requireNamespace("affxparser") || throw("Package not loaded: affxparser")
  affxparser::readPgfEnv(...)
}

.readPgfHeader <- function(...) {
  requireNamespace("affxparser") || throw("Package not loaded: affxparser")
  affxparser::readPgfHeader(...)
}

.readCcgHeader <- function(...) {
  requireNamespace("affxparser") || throw("Package not loaded: affxparser")
  affxparser::readCcgHeader(...)
}

.readCcg <- function(...) {
  requireNamespace("affxparser") || throw("Package not loaded: affxparser")
  affxparser::readCcg(...)
}

.writeCdf <- function(...) {
  requireNamespace("affxparser") || throw("Package not loaded: affxparser")
  affxparser::writeCdf(...)
}

.writeCdfHeader <- function(...) {
  requireNamespace("affxparser") || throw("Package not loaded: affxparser")
  affxparser::writeCdfHeader(...)
}

.writeCdfUnits <- function(...) {
  requireNamespace("affxparser") || throw("Package not loaded: affxparser")
  affxparser::writeCdfUnits(...)
}

.writeCdfQcUnits <- function(...) {
  requireNamespace("affxparser") || throw("Package not loaded: affxparser")
  affxparser::writeCdfQcUnits(...)
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Imports from aroma.light
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
.normalizeFragmentLength <- function(...) {
  requireNamespace("aroma.light") || throw("Package not loaded: aroma.light")
  aroma.light::normalizeFragmentLength(...)
}

.normalizeQuantile <- function(...) {
  requireNamespace("aroma.light") || throw("Package not loaded: aroma.light")
  aroma.light::normalizeQuantile(...)
}

.normalizeQuantileSpline <- function(...) {
  requireNamespace("aroma.light") || throw("Package not loaded: aroma.light")
  aroma.light::normalizeQuantileSpline(...)
}

.calibrateMultiscan <- function(...) {
  requireNamespace("aroma.light") || throw("Package not loaded: aroma.light")
  aroma.light::calibrateMultiscan(...)
}

.robustSmoothSpline <- function(...) {
  requireNamespace("aroma.light") || throw("Package not loaded: aroma.light")
  aroma.light::robustSmoothSpline(...)
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Imports for BiocGenerics
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
.annotation <- function(...) {
  requireNamespace("BiocGenerics") || throw("Package not loaded: BiocGenerics")
  BiocGenerics::annotation(...)
}

`.annotation<-` <- function(..., value) {
  requireNamespace("BiocGenerics") || throw("Package not loaded: BiocGenerics")
  ns <- getNamespace("BiocGenerics")
  get("annotation<-", mode="function", envir=ns)(..., value=value)
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Imports for Biobase
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
.assayData <- function(...) {
  requireNamespace("Biobase") || throw("Package not loaded: Biobase")
  Biobase::assayData(...)
}

`.assayData<-` <- function(..., value) {
  requireNamespace("Biobase") || throw("Package not loaded: Biobase")
  ns <- getNamespace("Biobase")
  get("assayData<-", mode="function", envir=ns)(..., value=value)
}

`.phenoData<-` <- function(..., value) {
  requireNamespace("Biobase") || throw("Package not loaded: Biobase")
  ns <- getNamespace("Biobase")
  get("phenoData<-", mode="function", envir=ns)(..., value=value)
}

`.protocolData<-` <- function(..., value) {
  requireNamespace("Biobase") || throw("Package not loaded: Biobase")
  ns <- getNamespace("Biobase")
  get("protocolData<-", mode="function", envir=ns)(..., value=value)
}

.featureNames <- function(...) {
  requireNamespace("Biobase") || throw("Package not loaded: Biobase")
  Biobase::featureNames(...)
}

`.featureNames<-` <- function(..., value) {
  requireNamespace("Biobase") || throw("Package not loaded: Biobase")
  ns <- getNamespace("Biobase")
  get("featureNames<-", mode="function", envir=ns)(..., value=value)
}

.sampleNames <- function(...) {
  requireNamespace("Biobase") || throw("Package not loaded: Biobase")
  Biobase::sampleNames(...)
}

`.sampleNames<-` <- function(..., value) {
  requireNamespace("Biobase") || throw("Package not loaded: Biobase")
  ns <- getNamespace("Biobase")
  get("sampleNames<-", mode="function", envir=ns)(..., value=value)
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Imports for oligo
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
.cleanPlatformName <- function(...) {
  requireNamespace("oligo") || throw("Package not loaded: oligo")
  oligo::cleanPlatformName(...)
}

.getPlatformDesign <- function(...) {
  requireNamespace("oligo") || throw("Package not loaded: oligo")
  oligo::getPlatformDesign(...)
}

.geometry <- function(...) {
  requireNamespace("oligo") || throw("Package not loaded: oligo")
  oligo::geometry(...)
}
