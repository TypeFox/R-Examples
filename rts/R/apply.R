# Author: Babak Naimi, naimi.b@gmail.com
# Date :  July 2012
# Version 1.0
# Licence GPL v3

if (!isGeneric("apply.daily")) {
  setGeneric("apply.daily", function(x, FUN, ...)
    standardGeneric("apply.daily"))
}

if (!isGeneric("apply.weekly")) {
  setGeneric("apply.weekly", function(x, FUN, ...)
    standardGeneric("apply.weekly"))
}

if (!isGeneric("apply.monthly")) {
  setGeneric("apply.monthly", function(x, FUN, ...)
    standardGeneric("apply.monthly"))
}

if (!isGeneric("apply.quarterly")) {
  setGeneric("apply.quarterly", function(x, FUN, ...)
    standardGeneric("apply.quarterly"))
}

if (!isGeneric("apply.yearly")) {
  setGeneric("apply.yearly", function(x, FUN, ...)
    standardGeneric("apply.yearly"))
}

setMethod("apply.daily", "RasterStackBrickTS",
          function(x, FUN, ...) {
            ep <- endpoints(x@time, "days")
            period.apply(x, ep, FUN, ...)
          })

setMethod("apply.weekly", "RasterStackBrickTS",
          function(x, FUN, ...) {
            ep <- endpoints(x@time, "weeks")
            period.apply(x, ep, FUN, ...)
          })


setMethod("apply.monthly", "RasterStackBrickTS",
          function(x, FUN, ...) {
            ep <- endpoints(x@time, "months")
            period.apply(x, ep, FUN, ...)
          })

setMethod("apply.quarterly", "RasterStackBrickTS",
          function(x, FUN, ...) {
            ep <- endpoints(x@time, "quarters")
            period.apply(x, ep, FUN, ...)
          })


setMethod("apply.yearly", "RasterStackBrickTS",
          function(x, FUN, ...) {
            ep <- endpoints(x@time, "years")
            period.apply(x, ep, FUN, ...)
          })

