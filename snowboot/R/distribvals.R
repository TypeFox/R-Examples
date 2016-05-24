distribvals <- function(distrib, vals, x) {
      if (any(vals == x)) {
            res <- distrib[which(vals == x)]
      } else {
            res <- 0
      }
      res
}
