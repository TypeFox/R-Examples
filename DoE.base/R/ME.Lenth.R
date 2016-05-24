ME.Lenth <- function(b, simulated=TRUE, alpha=NULL){
    s0 <- 1.5 * median(abs(b))
    cj <- as.numeric(b[abs(b) < 2.5 * s0])
    PSE <- 1.5 * median(abs(cj))
    m <- length(b)
    alphas <- as.numeric(colnames(crit.ME))
      if (is.null(alpha)) alpha <- alphas
      else{
         if (!is.numeric(alpha)) stop("alpha must be numeric")
         if (min(alpha)<0 || max(alpha)>1) stop("alpha must be between 0 and 1")
      }
      if (simulated){
        if (length(setdiff(alpha, alphas))>0 || m<7 || m>143){
        simulated <- FALSE
        message("simulated critical values not available for all requests, used conservative ones")
        }
      }
    ## cover simulated critical values, where available
    if (simulated){
      ME <- crit.ME[as.character(m),as.character(alpha)]*PSE
      SME <- crit.SME[as.character(m),as.character(alpha)]*PSE
    }
    else{
      gamma <- (1 + (1 - alpha)^(1/m))/2
      ME <- qt(1-alpha/2, m/3)*PSE
      SME <- qt(gamma, m/3)*PSE
      names(ME) <- names(SME) <- alpha
    }
    list(s0=s0, PSE=PSE, ME=ME, SME=SME)
  }

critsimul <- function(m, dfe, nsimul=10000, alpha=round(seq(0.25,0.01,-0.01),2), 
                      type="EM08", weight0=5){
## m the number of effects
## dfe the number of error df
## type "EM08" or "LW98"
## weight0 relevant for EM08 only
   ergs.CME <- rep(NA, m*dfe)
   ergs.CSME <- rep(NA, dfe)
   
   for (i in 1:nsimul){
      x <- rnorm(m)                   ## all effects inactive
      se <- sqrt(rchisq(1,dfe)/dfe)   ## estimated standard error
      if (type=="EM08") hilf <- CME.EM08(x, se, dfe, weight0=weight0, alpha=alpha, 
                    simulated=FALSE)
      else hilf <- CME.LW98(x, se, dfe, alpha=alpha, 
                    simulated=FALSE)
      hilf <- x/hilf$CPSE
      ergs.CME[((i-1)*m+1):(i*m)] <- hilf
      ergs.CSME[i] <- max(hilf)
   }
   CME.mult <- quantile(abs(ergs.CME), 1-alpha)
   CSME.mult <- quantile(abs(ergs.CSME), 1-alpha)
   names(CME.mult) <- alpha
   names(CSME.mult) <- alpha
   list(CME.mult = CME.mult, CSME.mult = CSME.mult)
}

CME.LW98 <- function(b, sterr, dfe, simulated=TRUE, alpha=NULL){
  ## sterr = sqrt(K*MSE), obtainable from linear model analysis
  ## dfe = residual df
  s0 <- 1.5 * median(abs(b))
  cj <- as.numeric(b[abs(b) < 2.5 * s0])
  m <- length(b)
  d <- m/3
  wPSE <- d / (d + dfe)
  CPSE <- sqrt((1.5 * median(abs(cj)))^2*wPSE + (1-wPSE)*sterr^2)
    alphas <- as.numeric(colnames(crit.ME))
      if (is.null(alpha)) alpha <- alphas
      else{
         if (!is.numeric(alpha)) stop("alpha must be numeric")
         if (min(alpha)<0 || max(alpha)>1) stop("alpha must be between 0 and 1")
      }
  if (simulated){
    hilf <- critsimul(m, dfe, type="LW98", alpha=alpha)
    CME <- hilf$CME.mult*CPSE
    CSME <- hilf$CSME.mult*CPSE
  }
  else{
    gamma <- (1 + (1 - alpha)^(1/m))/2
    CME <- qt(1-alpha/2, d+dfe)*CPSE
    CSME <- qt(gamma, d+dfe)*CPSE
    names(CME) <- names(CSME) <- alpha
  }
  list(s0=s0, CPSE=CPSE, CME=CME, CSME=CSME)
}

CME.EM08 <- function(b, sterr, dfe, simulated=TRUE, weight0=5, alpha=NULL){
  ## sterr = sqrt(K*MSE), obtainable from linear model analysis
  ## dfe = residual df
  ## weight0 <- multiplier for dfe in weighting sterr^2 for s0
  m <- length(b)
  d <- length(b)/3
  ws0 <- d / (d + weight0 * dfe)
  wPSE <- d / (d + dfe)
  s0 <- 1.5 * median(abs(b))
  s0 <- sqrt(ws0*s0^2 + (1-ws0)*sterr^2)
  cj <- as.numeric(b[abs(b) < 2.5 * s0])
  CPSE <- sqrt((1.5 * median(abs(cj)))^2*wPSE + (1-wPSE)*sterr^2)
    alphas <- as.numeric(colnames(crit.ME))
      if (is.null(alpha)) alpha <- alphas
      else{
         if (!is.numeric(alpha)) stop("alpha must be numeric")
         if (min(alpha)<0 || max(alpha)>1) stop("alpha must be between 0 and 1")
      }
  if (simulated){
    hilf <- critsimul(m, dfe, type="EM08", weight0=weight0, alpha=alpha)
    CME <- hilf$CME.mult*CPSE
    CSME <- hilf$CSME.mult*CPSE
  }
  else{
    gamma <- (1 + (1 - alpha)^(1/m))/2
    CME <- qt(1-alpha/2, d+dfe)*CPSE
    CSME <- qt(gamma, d+dfe)*CPSE
    names(CME) <- names(CSME) <- alpha
  }
  list(Cs0=s0, CPSE=CPSE, CME=CME, CSME=CSME)
}
