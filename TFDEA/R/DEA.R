#******************************************************************************
#
# Copyright Tom Shott, ETA, PSU, 2013, 2014, 2015
# Use granted under BSD license terms
#
# R TFDEA Package
#
# NOTE: See utility.R for description of naming nomenclature, use of .l & .b suffix
#
# DEA wraper function - public interface. Uppercase to distinguish from Benchmarking
#
#******************************************************************************
#
# Using functions from other files
#source("dea_internal.R")
#source("utility.R")

DEA <- function(x, y, rts="vrs", orientation="input",
                slack=TRUE, dual=FALSE,
                second="none", z=0,
                round=FALSE, debug=1){

  rts         <- .checkOption(rts,          "rts",          options.rts.l)
  orientation <- .checkOption(orientation,  "orientation",  options.orientation.l)

  slack       <- .checkOption(slack,        "slack",        TRUE)
  dual        <- .checkOption(dual,         "dual",         TRUE)

  second      <- .checkOption(second,       "second",       options.second.l)

  round       <- .checkOption(round,        "round",        TRUE)
  debug       <- .checkOption(debug,        "debug",        1)

  # Check that x & y are legal inputs & convert to standard dimension arrays
  x <- .checkData(x, "x")
  y <- .checkData(y, "y")

  # .checkDataGood now includes this check
  #   if (nrow(x) != nrow(y))
  #     stop("Number of DMU's in inputs != number of DMU's in outputs", call. = FALSE)
  .checkDataGood(x, y, debug = debug)

  # Check secondary optimization parms
  zn <- rep(0, nrow(x))
  second.b <- (second == "min" || second == "max")

  if (second.b){
    z <- .checkVector(z, "z")
    if (nrow(x) != length(z))
      stop("secondary data size (rows) must match number of DMU's", call. = FALSE)

    if (slack)
      stop("Can not use both second and slack option; set slack=FALSE", call. = FALSE)

    zn <- ifelse( rep(second == "min", nrow(x)),  -z,  z)
  }

  results <- .dea_internal(x, y, rts=rts, orientation=orientation,
                           slack=slack, dual=dual,
                           second.b=second.b, zn=zn,
                           round=round, debug=debug)

  return(results)
}
