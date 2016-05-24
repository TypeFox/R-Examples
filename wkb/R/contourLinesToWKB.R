contourLinesToWKB <- function(obj, endian = "little") {
  if(!inherits(obj, "list") || length(obj) < 1 ||
     !all(vapply(
       X = obj,
       FUN = function(subobj) {
         all(c("level", "x", "y") %in% names(subobj)) &&
           length(subobj[["level"]]) == 1 &&
           length(subobj[["x"]]) == length(subobj[["y"]])
       },
       FUN.VALUE = logical(1)
     ))
  ) {
   stop("obj must be a list of contour lines")
  }
  level <- vapply(X = obj, FUN = "[[", FUN.VALUE = numeric(1), "level")
  wkb <- XYListToWKBLineString(obj, endian)
  data.frame(level = level, contour = wkb)
}
