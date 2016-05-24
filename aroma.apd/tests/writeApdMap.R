library("aroma.apd")
library("R.utils") # Arguments

verbose <- Arguments$getVerbose(TRUE)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# 1. Scan for existing CEL files
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# a) Scan current directory for CEL files
files <- list.files(pattern="[.](cel|CEL)$")
files <- files[!file.info(files)$isdir]
if (length(files) > 0 && require("affxparser")) {
  # b) Corresponding APD filenames
  celNames <- files
  apdNames <- gsub(pattern, ".apd", files)


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # 1. Copy the probe intensities from a CEL to an APD file
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  for (kk in 1) {
    verbose && enter(verbose, "Reading CEL file #", kk)
    cel <- readCel(celNames[kk])
    verbose && exit(verbose)

    chipType <- cel$header$chiptype
    verbose && enter(verbose, "Getting read map for '", chipType, "'")
    mapType <- chipType
    mapFile <- paste(mapType, ".apm", sep="")
    if (file.exists(mapFile)) {
      verbose && enter(verbose, "Reading read map from APD map file")
      readMap <- readApdMap(mapFile)$map
      verbose && exit(verbose)
    } else {
      verbose && enter(verbose, "Generating read map from CDF file")
      cdfFile <- findCdf(chipType)
      readMap <- readCdfUnitsMap(cdfFile)
      writeApdMap(mapFile, map=readMap)
      verbose && exit(verbose)
    }
    verbose && exit(verbose)

    if (!file.exists(apdNames[kk])) {
      verbose && enter(verbose, "Calculating write map from read map")
      writeMap <- invertMap(readMap)
      verbose && exit(verbose)

      verbose && enter(verbose, "Creating APD file #", kk)
      writeApd(apdNames[kk], data=cel$intensities, chipType=chipType,
                                     mapType=mapType, writeMap=writeMap)
      verbose && exit(verbose)
      rm(writeMap)
    }

    verbose && enter(verbose, "Verifying APD file #", kk)
    apd <- readApd(apdNames[kk], readMap=readMap)
    verbose && exit(verbose)
    stopifnot(identical(apd$intensities, cel$intensities))

    rm(cel, mapType, mapFile, apd)
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # 2. Read a subset of the units
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  units <- c(1, 20:205)
  cel <- readCelUnits(celNames[1], units=units)
  apd <- readApdUnits(apdNames[1], units=units)
  stopifnot(identical(apd, cel))


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # 3. Benchmark reading with and without read map
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # To minimize the overhead from searching for and reading CDF
  # file when benchmarking, we read the CDF already here
  cdfFile <- findCdf(chipType)
  cdf <- readCdfCellIndices(cdfFile, units=units)

  cat("Benchmarks for read APD file:\n")
  t0 <- system.time({
    # Approximately 10 times faster than without a read map
    apd <- readApdUnits(apdNames[1], cdf=cdf, readMap=readMap)
  })[3]
  cat(sprintf("a) with read map by vector: %.3fs [ 1.00x]\n", t0))

  t <- system.time({
    apd <- readApdUnits(apdNames[1], cdf=cdf)
  })[3]
  cat(sprintf("b) with read map from file: %.3fs [%5.2fx]\n", t, t/t0))

  t <- system.time({
    # Force no read map to be used
    apd <- readApdUnits(apdNames[1], cdf=cdf, readMap=NULL)
  })[3]
  cat(sprintf("c) without read map       : %.3fs [%5.2fx]\n", t, t/t0))
} # if (length(files) > 0)
