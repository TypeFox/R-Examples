## ------------------------------------------------------------------ ##
## Authors: Luke Keele, Ohio State University
##          Jason Morgan, Ohio State University 
## ------------------------------------------------------------------ ##

## ------------------------------------------------------------------ ##
## Main function.
## Rt/Rc:   Response for treatment and control.
## Dt/Dc:   Dose received.
## Zt/Zc:   encouragement received (i.e., the instrumental variable).
## b:       Null hypothesis for beta (default = 0).
## ------------------------------------------------------------------ ##
iv_sens <- function(Rt, Rc, Dt, Dc, Gamma = 6, GammaInc = 1) {

  ## Length of the response vector for the treated.
  n <- length(Rt)
  Zt <- rep(1, length(Rt))
  Zc <- rep(0, length(Rt))
  
  ## Check that the response and dose vectors have the same length.
  if (length(Rc) != n | length(Dt) != n | length(Dc) != n |
      length(Zt) != n | length(Zc) != n) {
    stop("Length of response, dose, and encouragement vectors not equal.")
  }
  
  #Ensure test is in the same direction regardless of effect sign.
   est <- mean(Rt) - mean(Rc)
   
  if(est > 0){
    rt <- Rt
    rc <- Rc
      } else {
       rt <- Rc
       rc <- Rt
    }
  rm(Rt, Rc)
  #Rank sign test with IV adjusted outcomes
  b <- 0
  n <- length(rt)
  diff <- (Zt - Zc)*((rt - b*Dt) - (rc - b*Dc))
  diff <- diff[rt - rc != 0]
  S <- length(diff)
  ranks <- rank(abs(diff), ties.method="average")
  psi <- as.numeric(diff > 0)
  T <- sum(psi * ranks)
  gamma <- seq(1, Gamma, by=GammaInc)
  m <- length(gamma)
  pvals <- matrix(NA, m, 2)
  
  
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
 Instrument Due to Unobserved Factors \n"

  Obj <- list(Gamma = Gamma, GammaInc = GammaInc, pval = pval,
              msg = msg, bounds = bounds, note = note)
  class(Obj) <- c("rbounds", class(Obj))
  
  Obj
}   
 


