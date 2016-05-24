# ------------------------------------------------------------------------------
# Internal function 'surface.kernel'
#
# Author: Seong-Yun Hong <hong.seongyun@gmail.com>
# ------------------------------------------------------------------------------
surface.kernel <- function(coords, data, sigma, nrow, ncol, window, verbose) {
    
  if (verbose){
    begTime <- Sys.time(); fn <- match.call()[[1]]
    message(fn, ": kernel smoothing of the population data ...")
  }

  x <- coords[,1]; y <- coords[,2]
  
  for (i in 1:ncol(data)) {
    if (verbose)
      message(fn, ": processing column ", i)
    wgtXY <- cbind(rep(x, data[,i]), rep(y, data[,i]))
    tmp1 <- splancs::kernel2d(wgtXY, window, h0 = sigma, nx = ncol, ny = nrow, 
                              quiet = TRUE)
    
    # Transform the result to the format needed for spseg()
    tmp2 <- cbind(expand.grid(tmp1$x, tmp1$y), as.numeric(tmp1$z))
    colnames(tmp2) <- c("x", "y", "z")
    
    if (i == 1) {
      pixels <- as.matrix(tmp2[,1:2])
      values <- tmp2[,3]
    } else if (i > 1) {
      values <- cbind(values, tmp2[,3])
    }
  }

  # Remove points that are outside of the polygons
  outside <- which(is.na(values[,1]))
  if (length(outside) > 0) {
    pixels <- pixels[-outside,]
    values <- values[-outside,]    
  }

  if (verbose){
    tt <- as.numeric(difftime(Sys.time(), begTime, units = "sec"))
    message(fn, ": done! [", tt, " seconds]")
  }
  colnames(pixels) <- c("x", "y")
  colnames(values) <- colnames(data)
  list(coords = pixels, data = values)
}
