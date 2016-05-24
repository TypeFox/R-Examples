# Author: Babak Naimi, naimi.b@gmail.com
# Date :  July 2012
# Version 1.0
# Licence GPL v3

if (!isGeneric("endpoints")) {
  setGeneric("endpoints", function(x, on="months",k=1)
    standardGeneric("endpoints"))
}

setMethod("endpoints", "RasterStackBrickTS",
          function(x, on="months", k=1) {
            endpoints(x@time,on = on, k = k)
          }
          )

