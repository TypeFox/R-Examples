# ------------------------------------------------------------------------------
# Class 'SegDecomp'
#
# Author: Seong-Yun Hong <hong.seongyun@gmail.com>
# ------------------------------------------------------------------------------
setClass(Class = "SegDecomp", 
         slots = c(d = "numeric", coords = "matrix", data = "matrix",
                   proj4string = "CRS"))

setValidity(Class = "SegDecomp", 
  method = function(object) { 
    if (length(object@d) != 3)
      paste("length of 'd' must be three")
    else if (any(object@d < 0))
      paste("'d' must be positive")
    else if (sum(object@d) > 1)
      paste("sum of 'd' must be between 0 and 1")
    else if (ncol(object@coords) != 2 || !is.numeric(object@coords))
      paste("'coords' must be a numeric matrix with two columns, x and y")
    else if (ncol(object@data) < 2 || !is.numeric(object@data))
      paste("'data' must be a numeric matrix with at least two columns")
    else if (nrow(object@data) != nrow(object@coords))
      paste("'data' must have the same number of rows as 'coords'")
    else if (class(object@proj4string) != "CRS")
      paste("'proj4string' is not a valid CRS object")
    else
      TRUE
  })

SegDecomp <- function(d, coords, data, proj4string = CRS(as.character(NA))) {
  new("SegDecomp", d = d, coords = coords, data = data, 
      proj4string = proj4string)
}
