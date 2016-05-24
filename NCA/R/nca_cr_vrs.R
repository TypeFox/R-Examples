nca_cr_vrs <- 
function (loop.data, mpy, cutoff, bottleneck.x, fast.vrs=TRUE) {
  x <- loop.data$x
  y <- loop.data$y 
  
  # Transform the X and Y axes into positive values (for DEA)
  x.min <- min(x, 0)
  y.min <- min(y, 0)
  x <- x - x.min
  y <- y - y.min

  # Find the points on the ceiling ("PEERS")
  if (fast.vrs) {
    # Use out own implementation for finding peers
    peers <- p_vrs_peers(x, y)
  } else {
    # with optimal technical efficiency for DEA (vrs)
    vrs <- dea(x, y, RTS="vrs", ORIENTATION="graph")

    # Get the sorted, corrected peer matrix
    peers <- p_optimal_peers(vrs, x, y)
  }

  if (!is.vector(peers) && length(peers) > 2) {
    # Perform OLS through the peers
    x <- peers[,1] + x.min
    y <- peers[,2] + y.min
    line <- lm(y~x)
    
    intercept <- unname(coef(line)["(Intercept)"])
    slope     <- unname(coef(line)["x"])
    ceiling   <- p_ceiling(loop.data, slope, intercept)
    above     <- p_above(loop.data, slope, intercept)
  } else {
    ceiling   <- 0
    intercept <- NA
    slope     <- NA
    line      <- NULL
    above     <- 0
  }
  
  effect      <- ceiling / loop.data$scope
  ineffs      <- p_ineffs(loop.data, intercept, slope)
  bottleneck  <- p_bottleneck(loop.data, mpy, slope, intercept, cutoff, bottleneck.x)
  
  return(list(line=line,
              ceiling=ceiling, slope=slope, effect=effect,
              intercept=intercept, above=above, ineffs=ineffs,
              bottleneck=bottleneck))
}