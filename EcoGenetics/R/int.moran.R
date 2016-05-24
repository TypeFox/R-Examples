#' Moran internal.
#' 
#' @param Z Vector, matrix or data frame.
#' @param con Connection network.
#' @param nsim Number of Monte-Carlo simulations. 
#' @param alternative The alternative hypothesis. If "auto" is selected (default) the
#' program determines the hypothesis by difference between the median of the simulations
#' and the observed value. Other options are: "two.sided", "greater" and "less".
#' if test == cross, for the first interval (d == 0) the p and CI are computed with cor.test.
#' @param adjust.n Should be adjusted the number of individuals? (warning, this would
#' change variances)
#' @param plotit Should be generated a plot of the simulations?
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar}
#' 
#' @keywords internal

int.moran <- function(Z, con, nsim,
                      alternative, 
                      test = "permutation", adjust.n = FALSE, 
                      plotit) {
  
  N <- length(Z)
  wg <- int.check.con(con)
  
  #weight adjustment to number of connections
  if(adjust.n == TRUE) {
    colTRUE <- apply(wg, 1, sum)
    colTRUE[colTRUE != 0] <- 1
    Nc <- sum(colTRUE)
  } else  {
    Nc <- N
  }
  
  #Moran's I computation 
  z2 <- Z - mean(Z)
  SC <- drop(z2 %*% z2)
  VAR.Z <- SC / Nc
  W <- sum(wg)
  
  moranfun <- function(zc) {
    SCW <- wg %*% zc
    AUTOCOV.Z <- drop(zc %*% SCW) / W
    res <- AUTOCOV.Z / VAR.Z
    res
  }
  
  #observed value
  obs <- moranfun(z2)
  
  #Monte carlo replicates
  repsim <- numeric()
  for(i in 1:nsim) {
    samp <- sample(N)
    coefsup.mc <- z2[samp]
    repsim[i] <- moranfun(coefsup.mc)
  }
  
  
  #p value or CI computation
  random.m <- int.random.test(repsim = repsim, 
                              obs = obs,
                              nsim = nsim,
                              test = test, 
                              alternative = alternative)
  
  
  #listing results
  if(test == "permutation") {
    
    res <- list("analysis" = "Moran's I", 
                "alternative"= random.m$alter, 
                "observation" = round(random.m$obs, 4),
                "expectation" = round(random.m$exp, 4),
                "nsim" = nsim,
                "p.value" = round(random.m$p.val, 5),
                "quantile" = round(random.m$CI, 4))
    
  }else {
    
    res <- list("analysis" = "Moran's I",
                "observation" = round(random.m$obs, 4),
                "nsim" = nsim,
                "quantile" = round(random.m$CI, 4))
  }
  
  
  #plot
  if(plotit == TRUE) {
    hist(c(repsim, random.m$obs),
         xlab = "Moran's I", 
         main = "Monte Carlo test")
    abline(v = obs, col = "red")
    points(obs, 0, col = "green", 
           pch = 15, cex = 3.6)
    text(random.m$obs, 0, "obs")
  }
  
  res
  
} 
