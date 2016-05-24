ibgaout <- function(betaout, nsbp, alpha) {
  
  index <- 2 * (1 - alpha) ^(1 / (nsbp - 1)) - 1
  I_index <- qnorm((1 + index) / 2)
  
  over <- betaout[, 1] / betaout[, 2]
  cdfv <- pnorm(over)
  P_value <- max(min(1 - cdfv ^ (nsbp - 1)), min(1 - (1 - cdfv) ^ (nsbp - 1)))
  
  Powerv <- powerval(betaout, nsbp, test = "IBGA", alpha)
  
  list(index = I_index, pvalue = P_value, power = Powerv)
  
}

