p_ineffs <-
function (loop.data, intercept, slope) {
  if (is.na(slope) || slope <= 0) {
    return( list(x=NA, y=NA, abs=NA, rel=NA) )
  }

  x.c     <- (loop.data$y.high - intercept) / slope
  x.c.max <- min(loop.data$x.high, x.c)
  y.c.min <- slope * loop.data$x.low + intercept

  ineffs.x    <- (loop.data$x.high - x.c.max) / (loop.data$x.high - loop.data$x.low)
  ineffs.y    <- (y.c.min - loop.data$y.low)  / (loop.data$y.high - loop.data$y.low)
  ineffs.rel  <- ineffs.x + ineffs.y - ineffs.x * ineffs.y
  ineffs.abs  <- loop.data$scope * ineffs.rel

  return( list(x=ineffs.x * 100, y=unname(ineffs.y) * 100,
               abs=ineffs.abs, rel=ineffs.rel * 100) )
}

p_ineffs_ce <-
function (loop.data, peers) {
  # if there is only one peer, the ceiling zone is zero
  if (is.vector(peers) || length(peers) == 2) {
    return( list(x=NA, y=NA, abs=NA, rel=NA) )
  }

  x.c.max <- tail(peers, n=1)[1]
  y.c.min <- peers[1,2]

  ineffs.x    <- (loop.data$x.high - x.c.max) / (loop.data$x.high - loop.data$x.low)
  ineffs.y    <- (y.c.min - loop.data$y.low)  / (loop.data$y.high - loop.data$y.low)
  ineffs.rel  <- ineffs.x + ineffs.y - ineffs.x * ineffs.y
  ineffs.abs  <- loop.data$scope * ineffs.rel

  return( list(x=ineffs.x * 100, y=unname(ineffs.y) * 100,
               abs=ineffs.abs, rel=ineffs.rel * 100) )
}