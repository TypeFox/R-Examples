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

    if (!file.exists(apdNames[kk])) {
      verbose && enter(verbose, "Creating APD file #", kk)
      chipType <- cel$header$chiptype
      writeApd(apdNames[kk], data=cel$intensities, chipType=chipType)
      verbose && exit(verbose)
    }

    verbose && enter(verbose, "Verifying APD file #", kk)
    apd <- readApd(apdNames[kk])
    verbose && exit(verbose)
    stopifnot(identical(apd$intensities, cel$intensities))
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # 2. Read a subset of the units
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  units <- c(1, 20:205)
  cel <- readCelUnits(celNames[1], units=units)
  apd <- readApdUnits(apdNames[1], units=units)
  stopifnot(identical(apd, cel))


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # 3. The same, but stratified on PMs and MMs
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  apd <- readApdUnits(apdNames[1], units=units, stratifyBy="pmmm",
                                                addDimnames=TRUE)
} # if (length(files) > 0)
