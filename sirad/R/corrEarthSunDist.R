corrEarthSunDist <-
  function(i) {
    d <- 1 + 0.0334 * cos(0.01721*i - 0.0552)
    d
  }

