lambda.interp <- function(lambda, s) {
  ## interpolation of lambda according to s
  ## take lambda = sfrac * left + (1 - sfrac * right)
  
  if(length(lambda) == 1){ #degenerated case of one lambda given
    lens <- length(s)
    left <- rep(1, lens)
    right <- left
    sfrac <- rep(1, lens)
  } else {
    s[s > max(lambda)] <- max(lambda)
    s[s < min(lambda)] <- min(lambda)
    lenLambda <- length(lambda)
      sfrac <- (lambda[1]-s)/(lambda[1] - lambda[lenLambda])
      lambda <- (lambda[1] - lambda)/(lambda[1] - lambda[lenLambda])
      coord <- approx(lambda, seq(lambda), sfrac)$y
      left <- floor(coord)
      right <- ceiling(coord)
      sfrac <- (sfrac - lambda[right]) / (lambda[left] - lambda[right])
      sfrac[left==right] <- 1
    }
  list(left = left, right = right, frac = sfrac)
}
