powerval <- function(betaout, nsbp, test = c("IBGA", "LRT"), alpha) {
  
  effsize <- betaout[, 1] / betaout[, 2]
  
  if(test == "LRT") {
    simdata <- matrix(rnorm(nsbp * 10000, effsize, 1), nsbp, 10000)
    oversqneg <- (simdata ^ 2) * (simdata > 0)
    oversqpos <- (simdata ^ 2) * (simdata < 0)
    indneg <- (alphaf(nsbp, colSums(oversqneg)) < alpha) 
    indpos <- (alphaf(nsbp, colSums(oversqpos)) < alpha)
    Powerv <- mean(indneg * indpos)
  } else if(test == "IBGA") { 
    P_index <- qnorm((1 - alpha) ^ (1 / (nsbp - 1)))
    Powerv <- 1 - prod(pnorm(effsize + P_index)) - prod(pnorm(-effsize + P_index)) + 
      prod(pnorm(-effsize + P_index) - pnorm(-effsize - P_index))
    #reject <- rep(0, 1000)
    #for(i in 1 : 1000){
      #cdfv <- pnorm(simdata[, i])
      #P_value <- max(min(1 - cdfv ^ (nsbp - 1)), min(1 - (1 - cdfv) ^ (nsbp - 1)))
      #if(P_value < alpha) reject[i] <- 1
    #}
    #mean(reject)
  }
  
  Powerv
  
}

