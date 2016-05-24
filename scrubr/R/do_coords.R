# assign consitent lat/lon variables to be able to act on data inputs
# more easily - and save original column names to rename data on return

do_coords <- function(x, lat, lon) {
  #if (is.null(attr(x, "lat_var")) || is.null(attr(x, "lon_var"))) {
  x <- guess_latlon(x, lat, lon)
  if (is.null(attr(x, "lat_var_orig"))) attr(x, "lat_var_orig") <- lat
  if (is.null(attr(x, "lon_var_orig"))) attr(x, "lon_var_orig") <- lon
  # }
  return(x)
}
