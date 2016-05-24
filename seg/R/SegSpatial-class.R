# ------------------------------------------------------------------------------
# Class 'SegSpatial'
#
# Author: Seong-Yun Hong <hong.seongyun@gmail.com>
# ------------------------------------------------------------------------------
setClass(Class = "SegSpatial",
         slots = c(d = "numeric", r = "numeric", h = "numeric", p = "matrix"),
         contains = "SegLocal")

setValidity(Class = "SegSpatial", 
  method = function(object) {
    if (nrow(object@p) != ncol(object@p))
      paste("'p' must be a square matrix")
    else
      TRUE
  })

SegSpatial <- function(d, r, h, p, coords, data, env, 
  proj4string = CRS(as.character(NA))) {
  new("SegSpatial", d = d, r = r, h = h, p = p, coords = coords, data = data, 
      env = env, proj4string = proj4string)
}
