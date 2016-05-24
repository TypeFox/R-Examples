# 
# Helpers for grid building
# 
# 
# 
# Builds the grid, dispatching to the right function depending on provided
# grid type.
build_grid <- function(grid_type, coords, npts, pad, grid_opts) {
  
  grid.genfun.name <- paste0('build_grid_', grid_type)
  
  grid.genfun <- get(grid.genfun.name)
  return(as.data.frame(grid.genfun(coords, npts, pad, grid_opts)))
}

# Builds a set of equally-spaced points, taking pad into account. 
# This function is used internally by other grid-building functions.
build_grid_seed_onedim <- function(col, range, npts, pad) {
  seq(min(range[ ,col]) - pad, max(range[ ,col]) + pad, length.out = npts) 
}
