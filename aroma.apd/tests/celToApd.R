library("aroma.apd")
library("R.utils")  ## Arguments

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# 1. Scan for existing CEL files
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# a) Scan for CEL files
files <- list.files(pattern="[.](cel|CEL)$")
files <- files[!file.info(files)$isdir]
if (length(files) > 0 && require("affxparser")) {
  cat("Create an optimal read map for CEL file:", files[1], "\n")
  cdffile <- findCdf(readCelHeader(files[1])$chiptype)
  res <- cdfToApdMap(cdffile)

  cat("Converting CEL file to APD file:", files[1], "\n")
  apdfile <- celToApd(files[1])
  cat("Created APD file:", apdfile, "\n")
  file.remove(apdfile)

  cat("Converting CEL file to APD file with an optimized read map:", files[1], "\n")
  apdfile <- celToApd(files[1], mapType=res$mapType)
  cat("Created APD file:", apdfile, "\n")

  writeMap <- invertMap(res$readMap)
  for (file in files[-1]) {
    cat("Converting CEL file to APD file with an optimized read map:", file, "\n")
    apdfile <- celToApd(file, mapType=res$mapType, writeMap=writeMap)
    cat("Created APD file:", apdfile, "\n")
  }
} # if (length(files) > 0)
