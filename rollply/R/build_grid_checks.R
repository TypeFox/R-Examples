# 
# This function does the checks for the grid-builders

build_grid_check_vars <- function(coords, npts) { 
  
  # Check that all args are non-null
  if ( is.null(coords) || is.null(npts)) { 
    stop('One or more arguments is NULL')
  }
  
}