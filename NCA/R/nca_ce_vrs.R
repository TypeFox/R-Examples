nca_ce_vrs <-
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

  # if there is only one peer, the ceiling zone is zero
  ceiling <- 0
  if (!is.vector(peers) && length(peers) > 2) {
    # Transform the X and Y back
    peers[, 1] <- peers[, 1] + x.min
    peers[, 2] <- peers[, 2] + y.min

    for (i in 1:(nrow(peers)-1)) {
      part.a <- (peers[i+1,2] - peers[i,2]) * (peers[i+1,1] - peers[1,1])
      part.b <- 0.5 * (peers[i+1,2] - peers[i,2]) * (peers[i+1,1] - peers[i,1])
      ceiling <- ceiling + part.a - part.b
    }
  }
  effect      <- ceiling / loop.data$scope
  ineffs      <- p_ineffs_ce(loop.data, peers)
  bottleneck  <- p_bottleneck_vrs(loop.data, mpy, peers, cutoff, bottleneck.x)
  
  return(list(line=list(x + x.min, y + y.min),
              ceiling=ceiling, slope=NA, effect=effect,
              intercept=NA, above=0, ineffs=ineffs,
              bottleneck=bottleneck))
}