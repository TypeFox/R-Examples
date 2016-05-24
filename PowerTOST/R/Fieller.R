#------------------------------------------------------------------------------
# power / sample size of ratio of means with normality on original scale
# (based on Fieller's CI and Sasabuchi's test)
# Hauschke D., Kieser M., Diletti E. and Burke M.
# "Sample size determination for proving equivalence based on the ratio
#  of two means for normally distributed data"
# Stat. Med. 18(1) p93-105 (1999) 
# 
# Author: dlabes, started coding Jul 2010
#------------------------------------------------------------------------------

# the multivariate t- and normal distribution with arbitrary corr
# require(mvtnorm)

power.RatioF <- function(alpha=0.025, theta1=0.8, theta2, theta0=0.95, CV, CVb, 
                         n, design="2x2", setseed=TRUE)
{
  if (!(tolower(design) %in% c("2x2","parallel"))) {
    design <- "2x2"; warning("Design unknown. Will use 2x2 cross-over.\n",
        call.=FALSE, immediate.=TRUE)
  }
  if (missing(CV)) stop("CV must be given!")
  if (tolower(design)!="parallel"){
    if (missing(CVb)) stop("CVb(etween) must be given in case of cross-over!",
                           call.=FALSE)
  }
  if (missing(n))  stop("Number of subjects must be given!", call.=FALSE)
  if (n<=2) stop("Number of subjects must be >2!", call.=FALSE)
  
  if (missing(theta2)) theta2 <- 1/theta1
  
  # without set.seed() each run gives different power values
  # due to the implementation of pmvt() and pmvnorm(), 
  # especially if power is small
  if (setseed) set.seed(123456789)
  power <- .power.RatioF(alpha, theta1, theta2, theta0, CV, CVb, n, design)

  return(power)
}
# -------------------------------------------------------------------------
# power working horse without any input check
.power.RatioF <- function(alpha, theta1, theta2, theta0, CV, CVb, n, design)
{
  df <- n-2 #for both designs
  
  # non-centrality parms
  delta <- c(0,0)
  if (tolower(design)=="parallel"){
    delta[1] <- (theta0 - theta1)/CV/sqrt(2.*(1.0+theta1^2)/n)
    delta[2] <- (theta0 - theta2)/CV/sqrt(2.*(1.0+theta2^2)/n)
    # in case of theta0=theta1 or =theta2 and CV=0 NaN is produced
    # setting it to zero, is this correct?
    delta[is.nan(delta)] <- 0
    rho <- (1.0+theta1*theta2)/sqrt((1+theta1^2)*(1+theta2^2))
  } else {
    # cross-over
    CVw <- CV
    delta[1] <- (theta0 - theta1)/sqrt((CVw^2*(1.0+theta1^2) + 
                 CVb^2*(1.0-theta1)^2)/n)
    delta[2] <- (theta0 - theta2)/sqrt((CVw^2*(1.0+theta2^2) + 
                 CVb^2*(1.0-theta2)^2)/n)
    # in case of theta0=theta1 or =theta2 and CV=0 NaN is produced
    # setting it to zero, is this correct?
    delta[is.nan(delta)] <- 0
    
    rho <- (CVw^2*(1.0+theta1*theta2) + 
            CVb^2*(1.0+theta1*theta2-theta1-theta2)) /
            sqrt((CVw^2*(1.0+theta1^2) + CVb^2*(1.0-theta1)^2)*
                 (CVw^2*(1.0+theta2^2) + CVb^2*(1.0-theta2)^2))
  }
  # correlation
  corr <- diag(2)
  corr[1,2] <- corr[2,1] <- rho
  # upper limits of integration
  upper <- c(Inf,0)
  upper[2] <- -qt(1-alpha,df)

  # TODO: vectorized form of that
  if (df<4000){
    power <- pmvt(upper=upper, delta=delta, corr=corr, df=df, abseps=1e-9)
    upper[1] <- qt(1-alpha,df)
    power <- power - pmvt(upper=upper, delta=delta, corr=corr, df=df,
                          abseps=1e-9)
  } else {
    # large sample approx. multivariate normal distri.
    power <- pmvnorm(upper=upper, mean=delta, corr=corr)
    upper[1] <- qt(1-alpha,df)
    power <- power - pmvnorm(upper=upper, mean=delta, corr=corr)
  }
  
  # get rid of the attributes
  attributes(power) <- NULL
  # power<0 may happen due to numeric inaccuracies
  power <- ifelse(power<0, 0, power)
  return(power)
}
# ------------------------------------------------------------------------
# internal function: start of sample size search
.sampleN0.RatioF <- function(alpha=0.025, targetpower, theta1, theta2, theta0,
                             CV, CVb, design="2x2")
{
  # Hauschke D., Steinijans V. and Pigeot I.
  # "Bioequivalence studies in Drug Development"
  # Chapter 10.3.3
  # John Wiley & Sons, Chichester (2007)
  z1 <- qnorm(1-alpha)
  # value 0.04 corresponds roughly to log(0.96)
  # with lower values there can be many steps between 0.95 and 1
  if (abs(log(theta0))>0.04) z2 <- qnorm(targetpower) else {
    z2 <- qnorm(1-(1-targetpower)/2)
    theta0 <- 1
  }
  
  if (theta0<=theta1 | theta0>=theta2) stop("Ratio0 ",theta0,
      " must be between theta1/theta2 = ", theta1," / ",theta2," !")
  if (tolower(design)=="parallel"){
    n01 <- 2.0*(1+theta1^2)*(CV/(theta1-theta0))^2*(z1+z2)^2
    n02 <- 2.0*(1+theta2^2)*(CV/(theta2-theta0))^2*(z1+z2)^2
  } else {
  # cross-over
    CVw <- CV
    n01 <- (CVw^2*(1.0+theta1^2)+ CVb^2*(1.0-theta1)^2)*(z1+z2)^2 / 
           (theta1-theta0)^2
    n02 <- (CVw^2*(1.0+theta2^2)+CVb^2*(1.0-theta2)^2)*(z1+z2)^2 / 
           (theta2-theta0)^2
  } 
  n0 <- max(n01,n02)
  #make an even multiple of 2
  n0=ceiling(n0)
  n0 <- 2*trunc(n0/2)
  if (n0<4) n0 <- 4   # minimum sample size
  return(n0)
} # end function
# ---------------------------------------------------------------------------
# Sample size for a desired power of ratio of two means 
# with normality on original scale
# mainly intended for studies with clinical endpoints
# therefore 95% CIs had to be used and consequently alpha=0.025
sampleN.RatioF <- function(alpha=0.025, targetpower=0.8, theta1=0.8, theta2,
                           theta0=0.95, CV, CVb, design="2x2", print=TRUE, 
                           details=FALSE, imax=100, setseed=TRUE)
{
  if (!(tolower(design) %in% c("2x2","parallel"))) {
    design <- "2x2"
    warning("Design unknown. Will use 2x2 cross-over.",
            call.=FALSE, immediate.=TRUE)
  }
  if (missing(theta2)) theta2 <- 1/theta1
  if ( (theta0<=theta1) | (theta0>=theta2) ) {
    stop("Ratio ",theta0," not between margins ",theta1," / ",theta2,"!", 
         call.=FALSE)
  }
  
  if (missing(CV)) stop("Err: CV must be given!", call.=FALSE)
  
  if (tolower(design)!="parallel"){
    if (missing(CVb)) stop("CVb(etween) must be given in case of cross-over!",
                           call.=FALSE)
  } else CVb <- 0
  
  if (print) {
    cat("\n+++++++++++ Equivalence test - TOST +++++++++++\n")
    cat("    based on Fieller's confidence interval\n")
    cat("            Sample size estimation\n")
    cat("-----------------------------------------------\n")
    if (tolower(design=="parallel")) {
      cat("Study design: 2-group parallel\n")
    } else {
      cat("Study design: 2x2 cross-over\n")
    }
    cat("Ratio of means with normality on original scale\n")
    cat("alpha = ",alpha,", target power = ", targetpower,"\n", sep="")
    cat("BE margins        =",theta1,"...", theta2,"\n")
    if (tolower(design)=="parallel"){
      cat("Null (true) ratio = ",theta0,",  CV = ",CV,"\n", sep="")
    } else {
      cat("Null (true) ratio = ",theta0,",  CVw = ",CV,",  CVb = ",CVb,
          "\n", sep="")
    }
  }
  # Start
  n   <- .sampleN0.RatioF(alpha, targetpower, theta1, theta2, theta0, CV, CVb, 
                           design)
  if (setseed) set.seed(123456789)
  pow <- .power.RatioF(alpha, theta1, theta2, theta0, CV, CVb, n, design)
  if (details) {
    cat("\nSample size search\n")
    cat(" n     power\n")
    # do not print first too high
    if (pow<=targetpower) cat( n," ", formatC(pow, digits=6, format="f"),"\n")
  }
  
  steps <- 2
  iter  <- 0  # emergency brake
  # starting with too high power should be rare
  while (pow>targetpower) {
    if (n<=4) { # min number
      if (print & iter==0) cat( n," ", formatC(pow, digits=6, format="f"),"\n")
      break
    }
    n    <- n - steps
    iter <- iter+1
    if (setseed) set.seed(123456789)
    pow  <- .power.RatioF(alpha, theta1, theta2, theta0, CV, CVb, n, design)
    if (details) cat( n," ", formatC(pow, digits=6, format="f"),"\n")
    if (iter>imax) break 
  }
  while (pow<targetpower) {
    n    <- n+steps
    iter <- iter+1
    if (setseed) set.seed(123456789)
    pow <- .power.RatioF(alpha, theta1, theta2, theta0, CV, CVb, n, design)
    if (details) cat( n," ", formatC(pow, digits=6, format="f"),"\n")
    if (iter>imax) break 
  }
  if (print && !details) {
    cat("\nSample size\n")
    cat(" n     power\n")
    cat( n," ", formatC(pow, digits=6, format="f"),"\n")
  }
  #return
  res <- data.frame(design=design, alpha=alpha, CV=CV,CVb=CVb, theta0=theta0, 
                    theta1=theta1, theta2=theta2, n=n, power=pow, 
                    targetpower=targetpower)
  names(res) <-c("Design","alpha","CV","CVb","theta0","theta1","theta2",
                 "Sample size", "Achieved power", "Target power")
  
  if (print) return(invisible(res)) else return(res)
}
