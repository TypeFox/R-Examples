# ------------------------------------------------------------------------------
# Function 'spseg'
#
# Author: Seong-Yun Hong <hong.seongyun@gmail.com>
# ------------------------------------------------------------------------------
spseg <- function(x, data, method = "all", smoothing = "none", 
  nrow = 100, ncol = 100, window, sigma, useC = TRUE, negative.rm = FALSE, 
  tol = .Machine$double.eps, verbose = FALSE, ...) {

  # ----------------------------------------------------------------------------
  # STEP 1 Data preparation
  # ----------------------------------------------------------------------------
  if (verbose)
    tmp <- chksegdata(x, data)
  else
    tmp <- suppressMessages(chksegdata(x, data))
  coords <- tmp$coords; data <- tmp$data; proj4string <- tmp$proj4string
  smoothing <- 
    match.arg(smoothing, c("none", "kernel", "equal"), several.ok = FALSE)
  
  # ----------------------------------------------------------------------------
  # STEP 2 Estimate the data surface
  # ----------------------------------------------------------------------------
  if (smoothing == "equal") {
    tmp <- surface.equal(x, data, nrow, ncol, verbose)
    coords <- tmp$coords; data <- tmp$data
  }
  
  else if (smoothing == "kernel") {
    if (missing(window)) {
      x <- range(coords[,1]); y <- range(coords[,2])
      window <- matrix(c(x[1], y[1], 
                         x[1], y[2], 
                         x[2], y[2], 
                         x[2], y[1]), ncol = 2, byrow = TRUE)
    }
    
    if (missing(sigma))
      sigma <- min(bw.nrd(coords[,1]), bw.nrd(coords[,2]))
    
    tmp <- surface.kernel(coords, data, sigma, nrow, ncol, window, verbose)
    coords <- tmp$coords; data <- tmp$data
  }
    
  # ----------------------------------------------------------------------------
  # STEP 3 Calculate the population composition of each local environment
  # ----------------------------------------------------------------------------
  env <- localenv(coords, data, ...)
  env <- update(env, proj4string = CRS(proj4string))

  # ----------------------------------------------------------------------------
  # STEP 4 Compute the segregation indices
  # ----------------------------------------------------------------------------
  spatseg(env, method, useC, negative.rm, tol)
}
