#' Geary internal.
#' @param Z Vector, matrix or data frame.
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


int.geary <- function(Z, con, nsim,
                      alternative, test = "permutation", 
                      adjust.n = FALSE, 
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
  
  #Geary's C computation 
  wg <- as.numeric(wg)
  W <- sum(wg)
  z2 <- Z - mean(Z)
  SC <- sum(z2 %*% z2) * 2 * W
  VAR.Z <- SC / (Nc -1)
  mat <- as.matrix(dist(Z, upper = T))
  
  
  gearyfun <- function(distmat) {
    distmat <- as.numeric(distmat)
    num <- as.numeric(distmat ^ 2)
    res <- drop(wg %*% num) / VAR.Z
    res
  }
  
  #observed value
  obs <- gearyfun(mat)
  
  #Monte carlo replicates
  repsim <- numeric()
  for(i in 1:nsim) {
    samp <- sample(N)
    mat.mc <- mat[samp, samp]
    repsim[i] <- gearyfun(mat.mc)
  }
  
  
  #p value or CI computation
  random.m <- int.random.test(repsim = repsim, 
                              obs = obs,
                              nsim = nsim,
                              test = test, 
                              alternative = alternative)
  
  
  #listing results
  if(test == "permutation") {
    
    res <- list("analysis" = "Geary's I", 
                "alternative"= random.m$alter, 
                "observation" = round(random.m$obs, 4),
                "expectation" = round(random.m$exp, 4),
                "nsim" = nsim,
                "p.value" = round(random.m$p.val, 5),
                "quantile" = round(random.m$CI, 4)
    )
    
  }else {
    
    res <- list("analysis" = "Geary's C",
                "observation" = round(random.m$obs, 4),
                "nsim" = nsim,
                "quantile" = random.m$CI)
  }
  
  
  #plot
  if(plotit == TRUE) {
    hist(c(repsim, random.m$obs),
         xlab = "Geary's C", 
         main = "Monte Carlo test")
    abline(v = obs, col = "red")
    points(obs, 0, col = "green", 
           pch = 15, cex = 3.6)
    text(random.m$obs, 0, "obs")
  }
  
  res
  
}
