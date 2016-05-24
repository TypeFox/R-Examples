# calculating the extent ##works for Move and Movestack (both inherit SPDF)
setGeneric(".extcalc", function(obj, ext) standardGeneric(".extcalc"))
setMethod(f = ".extcalc", signature = c(obj = "SpatialPointsDataFrame", ext = "numeric"), 
          definition = function(obj, ext) {
            Range <- as.vector(c(obj@bbox[1, ], obj@bbox[2, ]))
            if (length(ext) == 1) ext <- rep(ext, 4)
            if (length(ext) == 2) ext <- rep(ext, each = 2)
            if (length(ext) == 4) {
              yRange <- c(Range[3] - abs(diff(Range[3:4]) * ext[3]), Range[4] + abs(diff(Range[3:4]) * 
                ext[4]))
              xRange <- c(Range[1] - abs(diff(Range[1:2]) * ext[1]), Range[2] + abs(diff(Range[1:2]) * 
                ext[2]))
            } else {
              stop("The ext argument must be a vector of 1, 2 or 4 numbers")
            }
            return(c(xRange, yRange))
          })
