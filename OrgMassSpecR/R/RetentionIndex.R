RetentionIndex <- function(n, target, preceding, following) {

  retention.index <- 100 * n + 100 * (target - preceding) / (following - preceding)
  return(retention.index)

}

