library("aroma.apd")

# Float precision
.Machine$float.eps <- (2^((8-4)*8)*.Machine$double.eps)
tol <- .Machine$float.eps ^ 0.5


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# 1. Create an Affymetrix Probe Signal (APD) file for a
#    'Mapping50K_Hind240' with 1600-by-1600 probes.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
chipType <- "Mapping50K_Hind240"
nbrOfCells <- 1600^2

pathname <- paste(tempfile(), "apd", sep=".")
createApd(pathname, nbrOfCells=nbrOfCells, chipType=chipType)

# File size
cat("File name:", pathname, "\n")
cat("File size:", file.info(pathname)$size, "bytes\n")
cat("APD header:\n")
header <- readApdHeader(pathname)
print(header)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# 2. Update the signals for a subset of probes
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cells <- c(1, 1100:1120, 220:201, 998300:999302)
signals <- log(cells+1, base=2)  # Fake signals
updateApd(pathname, indices=cells, data=signals)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# 3. Re-read the signals for a subset of probes
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
apd <- readApd(pathname, indices=cells)

# Signals in APD files are stored as floats (since this is
# the precision in CEL files).
stopifnot(all.equal(signals, apd$intensities, tolerance=tol))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# 4. Re-read the signals for a subset of probes
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if (require("affxparser")) {
  cdfFile <- findCdf(chipType)
  if (length(cdfFile) > 0) {

    apd <- readApdUnits(pathname, units=100:104)

    # Sample new data (with one decimal precision)
    apd2 <- lapply(apd, function(unit) {
      lapply(unit, function(groups) {
        n <- length(groups$intensities)
        values <- as.integer(runif(n, max=655350))/10
        list(intensities=values)
      })
    })

    # Update APD file with new data
    updateApdUnits(pathname, units=100:104, data=apd2)

    # Re-read data to verify correctness
    apd <- readApdUnits(pathname, units=100:104)

    # Signals in APD files are stored as floats (since this is
    # the precision in CEL files).
    stopifnot(all.equal(apd2, apd, tolerance=tol))
  } # if (length(cdfFile) > 0 ...)
} # if (require("affxparser"))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# 4. Clean up
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
file.remove(pathname)
