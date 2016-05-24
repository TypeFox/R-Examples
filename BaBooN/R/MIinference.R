# Inference via Rubin's Combining Rules 
# Version:       0.1-6
# Date:     2011-02-24
# Author:         F.M.

###################################################################
 
MI.inference <- function(thetahat, varhat.thetahat, alpha = 0.05)
{
    M <- length(thetahat)
    if (length(varhat.thetahat) != M) {
      stop("Different length for 'thetahat' and 'varhat.thetahat'!\n")}
    lambda <- 1 - (alpha/2)
    MIestimate <- mean(thetahat)
    B <- var(thetahat)  
    W <- mean(varhat.thetahat)
    total <- W + (1 + 1/M) * B
    DF <- (M - 1) * (1 + W/((1 + 1/M) * B))^2
    CI.low <- MIestimate - qt(lambda, DF) * sqrt(total)
    CI.up <- MIestimate + qt(lambda, DF) * sqrt(total)
    x <- list(MI.Est = MIestimate, MI.Var = total, CI.low = CI.low,
              CI.up = CI.up, BVar = B, WVar = W)
    return(x)
}


###################################################################
