psens <- function(x, y = NULL, Gamma = 6, GammaInc = 1){
  
  if (class(x)!="Match") {
    trt <- x
    ctrl <- y
  }
  else if(x$est > 0){
    ctrl <-x$mdata$Y[x$mdata$Tr==0]
    trt <- x$mdata$Y[x$mdata$Tr==1]
  } else {
    ctrl <- x$mdata$Y[x$mdata$Tr==1]
    trt <- x$mdata$Y[x$mdata$Tr==0]	
  }
  
  gamma <- seq(1, Gamma, by=GammaInc)
  m <- length(gamma)
  pvals <- matrix(NA, m, 2)
  diff <- trt - ctrl
  S <- length(diff)
  diff <- diff[diff != 0]
  ranks <- rank(abs(diff), ties.method="average")
  psi <- as.numeric(diff > 0)
  T <- sum(psi * ranks)
  
      for(i in 1:m) {
      p.plus <- gamma[i]/(1 + gamma[i])
      p.minus <- 1/(1+gamma[i])
      E.T.plus <- sum(ranks*p.plus)
      V.T <- sum(ranks^2 * p.plus*(1-p.plus))
      E.T.minus <- sum(ranks*p.minus)
      
      ## Normal Approximation
      z.plus <- (T - E.T.plus)/sqrt(V.T)
      z.minus <- (T - E.T.minus)/sqrt(V.T)
      p.val.up <- 1 - pnorm(z.plus)
      p.val.low <- 1 - pnorm(z.minus)
      pvals[i,1] <- round(p.val.low, digits=4)
      pvals[i,2] <- round(p.val.up, digits=4)
    }
  
  pval <- pvals[1,1]
  bounds <- data.frame(gamma, pvals)
  names(bounds) <- c("Gamma", "Lower bound", "Upper bound")

  msg <- "Rosenbaum Sensitivity Test for Wilcoxon Signed Rank P-Value \n"
  note <- "Note: Gamma is Odds of Differential Assignment To
 Treatment Due to Unobserved Factors \n"

  Obj <- list(Gamma = Gamma, GammaInc = GammaInc, pval = pval,
              msg = msg, bounds = bounds, note = note)
  class(Obj) <- c("rbounds", class(Obj))
  
  Obj
}
