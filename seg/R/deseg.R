# ------------------------------------------------------------------------------
# Function 'deseg'
#
# Author: Seong-Yun Hong <hong.seongyun@gmail.com>
# ------------------------------------------------------------------------------
deseg <- function(x, data, smoothing = "kernel", nrow = 100, ncol = 100, 
  window, sigma, verbose = FALSE) {

  # ----------------------------------------------------------------------------
  # STEP 1 Data preparation
  # ----------------------------------------------------------------------------
  if (verbose)
    tmp <- chksegdata(x, data)    
  else
    tmp <- suppressMessages(chksegdata(x, data))
  coords <- tmp$coords; data <- tmp$data; proj4string <- tmp$proj4string
  
  smoothing <- 
    match.arg(smoothing, "kernel", several.ok = FALSE)
      
  # ----------------------------------------------------------------------------
  # STEP 2 Estimate the data surface
  # ----------------------------------------------------------------------------
  if (smoothing == "kernel") {
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
  }
    
  # ----------------------------------------------------------------------------
  # STEP 3 Calculate the index
  # ----------------------------------------------------------------------------
  v <- .decomp(tmp$data)
  SL <- v - .decompL(tmp$data)
  SQ <- 1 / ncol(tmp$data)
  SC <- v - (SL + SQ)
  
  SegDecomp(d = c(SL, SC, SQ), tmp$coords, tmp$data, CRS(proj4string))
}
