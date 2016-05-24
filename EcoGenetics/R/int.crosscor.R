
#' Cross correlation. Internal.
#' @param Z vector, matrix or data frame.
#' with variable/s (in matrix or data frame formats, variables in columns).
#' @param Y Vector with the second variable for cross-correlograms construction.
#' if Z has multiple variables, the program will compute the cross correlograms 
#' between each and Y.
#' @param con Connection network.
#' @param nsim Number of Monte-Carlo simulations. 
#' @param alternative The alternative hypothesis. If "auto" is selected (default) the
#' program determines the hypothesis by difference between the median of the simulations
#' and the observed value. Other options are: "two.sided", "greater" and "less".
#' if test == cross, for the first interval (d== 0) the p and CI are computed with cor.test.
#' @param adjust.n Should be adjusted the number of individuals? (warning, this would
#' change variances)
#' @param plotit Should be generated a plot of the simulations?
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar}
#' @keywords internal


int.crosscor <- function(Z, Y, con, nsim,
                         alternative, test = "permutation", 
                         adjust.n = FALSE, 
                         plotit) {
  
  
  N <- length(Z)
  
  wg <- int.check.con(con)
  if(!isSymmetric(wg)) {
    stop("Non symmetric weight matrix. Weights must be symmetric.")
  }
  
  if(adjust.n == TRUE) {
    colTRUE <- apply(wg, 1, sum)
    colTRUE[colTRUE != 0] <- 1
    Nc <- sum(colTRUE)
  } else  {
    Nc <- N
  }
  
  
  z2 <- Z - mean(Z)
  y2 <- Y - mean(Y)
  W <- sum(wg)
  SCXX <- drop(z2 %*% z2)
  SCYY <- drop(y2 %*% y2)
  VAR.XY <- sqrt(SCXX * SCYY) / Nc
  Mi <- outer(z2, y2)
  
  crossfun  <- function(M) {
    AUTOCOV.XY <- sum(wg * M) / W
    res <- AUTOCOV.XY / VAR.XY
    res
  }
  
  obs <- crossfun(Mi)
  
  repsim <- numeric()
  
  for(i in 1:nsim) {
    samp <- sample(N)
    repsim[i] <- crossfun(Mi[samp, samp])  #permutation is conditioned to the each dependent pair (z2,y2)
  }					 		                           #symmetric weights required!
  random.m <- int.random.test(repsim = repsim, obs = obs, 
                              test = test,
                              nsim = nsim,
                              alternative = alternative)
  
  #monte carlo
  
  if(test == "permutation") {
    
    res <- list("analysis" = "bivariate Moran's I", 
                "alternative"=   random.m$alter, 
                "observation" = round(random.m$obs, 4),
                "expectation" = round(random.m$exp, 4),
                "nsim" = nsim, "p.value" = round(random.m$p.val, 5),
                "quantile" = round(random.m$CI, 4))
    
  }else {
    
    res <- list("analysis" = "Bivariate Moran's I", 
                "observation" = round(random.m$obs, 4),
                "nsim" = nsim,
                "quantile" = random.m$CI)
  }
  
  
  if(plotit == TRUE) {
    hist(c(repsim, random.m$obs), 
         xlab = "Bivariate Moran's Ixy", 
         main = "Monte Carlo test")
    abline(v = obs, col = "red")
    points(obs, 0, col = "green", pch = 15, cex = 3.6)
    text(random.m$obs, 0, "obs")
  }
  
  res
}
