reassign <- function(x) {
  if (!is.null(attr(x, "lat_var_orig"))) {
    x <- rename(x, c('latitude' = attr(x, "lat_var_orig")))
    x <- rename(x, c('longitude' = attr(x, "lon_var_orig")))
  }
  return(x)
}
