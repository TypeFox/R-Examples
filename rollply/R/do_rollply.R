# Do the actual job of computing the result
# 

do_rollply <- function(current_coords,   # the coordinates need to be the first arg
                       coords,
                       wdw.size,         # wdw size
                       .data,
                       fun,
                       ...) {
  
  # Get subset of values required
  lookup <- lookup_in_wdw(as.matrix(current_coords), coords, wdw.size)
  
  # Return function applied to that subset
  fun(.data[lookup, ],...)
}
