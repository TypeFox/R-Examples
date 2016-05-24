# Fits a circular von Mises kernel density
# Just a wrapper for the C function densRad
# NO INPUT CHECKS : USE WITH CARE!

# x is a sample of observations, time in radians
# grid is locations in (0, 2*pi] where the density will be calculated
# bw is the "band width" = concentration of the von Mises kernel
# Returns: a vector of densities estimated at the locations in 'grid'.

densityFit <-
function(x, grid, bw) {
  n <- length(x)
  nxpts <- length(grid)
  dens <- .C("densRad", as.double(x), as.integer(n),
    as.double(grid), as.integer(nxpts), as.double(bw),
    result = double(length(grid)))
  dens[['result']]
}
