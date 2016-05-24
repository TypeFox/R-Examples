MSClaio2008 <- function(sample, dist=c("NORM","LN","GUMBEL","EV2","GEV","P3","LP3"), crit=c("AIC", "AICc", "BIC", "ADC")) {
 # sample = sample
 # dist = distribution ("NORM","LN","EV1","EV2","GEV","GAM","LP3")
 # crit = criterion ("AIC", "AICc", "BIC", "ADC")
 MSC <- NULL
 MSC$sample <- sample
 MSC$dist <- dist
 MSC$crit <- crit
 n <- length(sample)
 nd <- length(dist)
 if (any("AIC"==crit)) {
  MSC$AIC <- rep(NA, nd)
  for (i in 1:nd) {
   MSC$AIC[i] <- .AIC(sample, dist[i])
  }
  names(MSC$AIC) <- dist
  MSC$AICdist <- dist[which.min(MSC$AIC)]
 }
 if (any("AICc"==crit)) {
  MSC$AICc <- rep(NA, nd)
  for (i in 1:nd) {
   MSC$AICc[i] <- .AICc(sample, dist[i])
  }
  names(MSC$AICc) <- dist
  MSC$AICcdist <- dist[which.min(MSC$AICc)]
 }
 if (any("BIC"==crit)) {
  MSC$BIC <- rep(NA, nd)
  for (i in 1:nd) {
   MSC$BIC[i] <- .BIC(sample, dist[i])
  }
  names(MSC$BIC) <- dist
  MSC$BICdist <- dist[which.min(MSC$BIC)]
 }
 if (any("ADC"==crit)) {
  MSC$ADC <- rep(NA, nd)
  for (i in 1:nd) {
   MSC$ADC[i] <- .ADC(sample, dist[i])
  }
  names(MSC$ADC) <- dist
  MSC$ADCdist <- dist[which.min(MSC$ADC)]
 }
 MSC$call <- match.call()
 class(MSC) <- "MSClaio2008"
 return(MSC)
}


#   MSClaio2008 <- function(sample, dist=c("NORM","LN","GUMBEL","EV2","GEV","P3","LP3"), crit=c("AIC", "AICc", "BIC", "ADC")) {
#    # sample = sample
#    # dist = distribution ("NORM","LN","EV1","EV2","GEV","GAM","LP3")
#    # crit = criterion ("AIC", "AICc", "BIC", "ADC")
#    MSC <- NULL
#    MSC$sample <- sample
#    MSC$dist <- dist
#    MSC$crit <- crit
#    n <- length(sample)
#    nd <- length(dist)
#    if (any("AIC"==crit)||any("AICc"==crit)||any("BIC"==crit)) {
#     lnMLs <- rep(NA, nd)
#     ps <- rep(2, nd)
#     for (i in 1:nd) {
#      lnMLs[i] <- .lnML(sample, dist[i])
#      if (any(c("GEV","GAM","P3","LP3") == dist[i])) ps[i] <- 3
#     }
#    }
#    if (any("AIC"==crit)) {
#     MSC$AIC <- -2*lnMLs + 2*ps
#     names(MSC$AIC) <- dist
#     MSC$AICdist <- dist[which.min(MSC$AIC)]
#     #MSC$AICpar <- ML_estimation(sample, MSC$AICdist)
#     #names(MSC$AICpar) <- paste("par", c(1:ps[which.min(MSC$AIC)]), sep="")
#    }
#    if (any("AICc"==crit)) {
#     MSC$AICc <- -2*lnMLs + 2*ps*(n/(n-ps-1))
#     names(MSC$AICc) <- dist
#     MSC$AICcdist <- dist[which.min(MSC$AICc)]
#     #MSC$AICcpar <- ML_estimation(sample, MSC$AICcdist)
#     #names(MSC$AICcpar) <- paste("par", c(1:ps[which.min(MSC$AICc)]), sep="")
#    }
#    if (any("BIC"==crit)) {
#     MSC$BIC <- -2*lnMLs + log(n)*ps
#     names(MSC$BIC) <- dist
#     MSC$BICdist <- dist[which.min(MSC$BIC)]
#     #MSC$BICpar <- ML_estimation(sample, MSC$BICdist)
#     #names(MSC$BICpar) <- paste("par", c(1:ps[which.min(MSC$BIC)]), sep="")
#    }
#    if (any("ADC"==crit)) {
#     MSC$ADC <- rep(NA, nd)
#     for (i in 1:nd) {
#      MSC$ADC[i] <- .ADC(sample, dist[i])
#     }
#     names(MSC$ADC) <- dist
#     MSC$ADCdist <- dist[which.min(MSC$ADC)]
#     #MSC$ADCpar <- ML_estimation(sample, MSC$ADCdist)
#     #names(MSC$ADCpar) <- paste("par", c(1:ps[which.min(MSC$ADC)]), sep="")
#    }
#    MSC$call <- match.call()
#    class(MSC) <- "MSClaio2008"
#    return(MSC)
#   }


# ----------------------------------------------------------- #

#   print.MSClaio2008 <- function(x, digits=max(3, getOption("digits") - 3), ...) {
#    cat("\nCall:\n", deparse(x$call), sep = "")
#       if (length(x$AIC)) {
#           cat("\n\n------------------------\nAkaike Information Criterion:\n\n")
#           print.default(format(x$AIC, digits=digits), print.gap=2, quote=FALSE)
#           cat("\n  Selected distribution: ", deparse(x$AICdist), "\n")
#           #print.default(format(x$AICpar, digits=digits), print.gap=2, quote=FALSE)
#       }
#       if (length(x$AICc)) {
#           cat("\n\n------------------------\nCorrected Akaike Information Criterion:\n\n")
#           print.default(format(x$AICc, digits=digits), print.gap=2, quote=FALSE)
#           cat("\n  Selected distribution: ", deparse(x$AICcdist), "\n")
#           #print.default(format(x$AICcpar, digits=digits), print.gap=2, quote=FALSE)
#       }
#       if (length(x$BIC)) {
#           cat("\n\n------------------------\nBayesian Information Criterion:\n\n")
#           print.default(format(x$BIC, digits=digits), print.gap=2, quote=FALSE)
#           cat("\n  Selected distribution: ", deparse(x$BICdist), "\n")
#           #print.default(format(x$BICpar, digits=digits), print.gap=2, quote=FALSE)
#       }
#       if (length(x$ADC)) {
#           cat("\n\n------------------------\nAnderson-Darling Criterion:\n\n")
#           print.default(format(x$ADC, digits=digits), print.gap=2, quote=FALSE)
#           cat("\n  Selected distribution: ", deparse(x$ADCdist), "\n")
#           #print.default(format(x$ADCpar, digits=digits), print.gap=2, quote=FALSE)
#       }
#       cat("\n")
#       invisible(x)
#   }

print.MSClaio2008 <- function(x, digits=max(3, getOption("digits") - 3), ...) {
 cat("\n")
 if (length(x$AIC)) {
  cat("------------------------\nAkaike Information Criterion (AIC):\n")
  print.default(format(x$AIC, digits=digits), print.gap=2, quote=FALSE)
 }
 if (length(x$AICc)) {
  cat("------------------------\nCorrected Akaike Information Criterion (AICc):\n")
  print.default(format(x$AICc, digits=digits), print.gap=2, quote=FALSE)
 }
 if (length(x$BIC)) {
  cat("------------------------\nBayesian Information Criterion (BIC):\n")
  print.default(format(x$BIC, digits=digits), print.gap=2, quote=FALSE)
 }
 if (length(x$ADC)) {
  cat("------------------------\nAnderson-Darling Criterion (ADC):\n")
  print.default(format(x$ADC, digits=digits), print.gap=2, quote=FALSE)
 }
 #cat("------------------------\nCall:\n", deparse(x$call), sep = "")
 cat("\n")
 invisible(x)
}



summary.MSClaio2008 <- function(object, ...) {
 x <- object
 cat("\nTested distributions:\n")
 print.default(x$dist, print.gap=2, quote=FALSE)
 cat("------------------------\n")

 Sdist <- NULL
 cont=1
 if (length(x$AICdist)) {
  Sdist <- c(Sdist, x$AICdist)
  names(Sdist)[cont] <- "AIC"
  cont <- cont+1
 }
 if (length(x$AICcdist)) {
  Sdist <- c(Sdist, x$AICcdist)
  names(Sdist)[cont] <- "AICc"
  cont <- cont+1
 }
 if (length(x$BICdist)) {
  Sdist <- c(Sdist, x$BICdist)
  names(Sdist)[cont] <- "BIC"
  cont <- cont+1
 }
 if (length(x$ADCdist)) {
  Sdist <- c(Sdist, x$ADCdist)
  names(Sdist)[cont] <- "ADC"
  cont <- cont+1
 }
 cat("Chosen distributions:\n")
 print.default(Sdist, print.gap=2, quote=FALSE)
 cat("whose Maximum-Likelihood parameters are:\n")

 Sdist <- unique(Sdist)
 if (any(Sdist=="NORM")) {
  T <- ML_estimation(x$sample, dist="NORM")
  cat("NORM parameters of x:", T, "\n", sep="  ")
 }
 if (any(Sdist=="LN")) {
  T <- ML_estimation(log(x$sample), dist="NORM")
  cat("NORM parameters of log(x):", T, "\n", sep="  ")
 }
 if (any(Sdist=="EV1")||any(Sdist=="GUMBEL")) {
  T <- ML_estimation(x$sample, dist="EV1")
  cat("GUMBEL parameters of x:", T, "\n", sep="  ")
 }
 if (any(Sdist=="EV2")) {
  T <- ML_estimation(log(x$sample), dist="EV1")
  cat("GUMBEL parameters of log(x):", T, "\n", sep="  ")
 }
 if (any(Sdist=="GEV")) {
  T <- ML_estimation(x$sample, dist="GEV")
  cat("GEV parameters of x:", T, "\n", sep="  ")
 }
 if (any(Sdist=="P3")||any(Sdist=="GAM")) {
  T <- ML_estimation(x$sample, dist="P3")
  cat("P3 parameters of x:", T, "\n", sep="  ")
 }
 if (any(Sdist=="LP3")) {
  T <- ML_estimation(log(x$sample), dist="P3")
  cat("P3 parameters of log(x):", T, "\n", sep="  ")
 }
 cat("\n")
}


# ----------------------------------------------------------- #

plot.MSClaio2008 <- function(x, ...) {
 lognormplot(x$sample, line=FALSE, col="blue", ...) 
 Sdist <- NULL
 if (length(x$AICdist)) Sdist <- c(Sdist, x$AICdist)
 if (length(x$AICcdist)) Sdist <- c(Sdist, x$AICcdist)
 if (length(x$BICdist)) Sdist <- c(Sdist, x$BICdist)
 if (length(x$ADCdist)) Sdist <- c(Sdist, x$ADCdist)
 Sdist <- unique(Sdist)
 cont=1
 lableg=NULL
 colleg=NULL
 if (any(x$dist=="NORM")) {
  T <- ML_estimation(x$sample, dist="NORM")
  q <- qnorm(seq(.001,.999,by=.001),mean=T[1],sd=T[2])
  if (any(Sdist=="NORM")) colore=1 else colore="gray"
  normpoints(q, type="l", lty=cont, col=colore)
  cont <- cont+1
  lableg <- c(lableg, "NORM")
  colleg <- c(colleg, colore)
 }
 if (any(x$dist=="LN")) {
  T <- ML_estimation(log(x$sample), dist="NORM")
  q <- exp(qnorm(seq(.001,.999,by=.001),mean=T[1],sd=T[2]))
  if (any(Sdist=="LN")) colore=1 else colore="gray"
  normpoints(q, type="l", lty=cont, col=colore)
  cont <- cont+1
  lableg <- c(lableg, "LN")
  colleg <- c(colleg, colore)
 }
 if (any(x$dist=="EV1")||any(x$dist=="GUMBEL")) {
  T <- ML_estimation(x$sample, dist="EV1")
  q <- T[1] - T[2]*log(-log(seq(.001,.999,by=.001)))
  if (any(Sdist=="EV1")||any(Sdist=="GUMBEL")) colore=1 else colore="gray"
  normpoints(q, type="l", lty=cont, col=colore)
  cont <- cont+1
  lableg <- c(lableg, "GUMBEL")
  colleg <- c(colleg, colore)
 }
 if (any(x$dist=="EV2")) {
  T <- ML_estimation(log(x$sample), dist="EV1")
  q <- exp(T[1] - T[2]*log(-log(seq(.001,.999,by=.001))))
  if (any(Sdist=="EV2")) colore=1 else colore="gray"
  normpoints(q, type="l", lty=cont, col=colore)
  cont <- cont+1
  lableg <- c(lableg, "EV2")
  colleg <- c(colleg, colore)
 }
 if (any(x$dist=="GEV")) {
  T <- ML_estimation(x$sample, dist="GEV")
  #q <- T[1] + (T[2]/T[3])*(1 - (-log(seq(.001,.999,by=.001)))^T[3])
  q <- invF.GEV(seq(.001,.999,by=.001), T[1], T[2], T[3])
  if (any(Sdist=="GEV")) colore=1 else colore="gray"
  normpoints(q, type="l", lty=cont, col=colore)
  cont <- cont+1
  lableg <- c(lableg, "GEV")
  colleg <- c(colleg, colore)
 }
 if (any(x$dist=="P3")||any(x$dist=="GAM")) {
  T <- ML_estimation(x$sample, dist="P3")
  q <- invF.gamma(seq(.001,.999,by=.001), T[1], T[2], T[3])
  if (any(Sdist=="P3")||any(Sdist=="GAM")) colore=1 else colore="gray"
  normpoints(q, type="l", lty=cont, col=colore)
  cont <- cont+1
  lableg <- c(lableg, "P3")
  colleg <- c(colleg, colore)
 }
 if (any(x$dist=="LP3")) {
  T <- ML_estimation(log(x$sample), dist="P3")
  q <- exp(invF.gamma(seq(.001,.999,by=.001), T[1], T[2], T[3]))
  if (any(Sdist=="LP3")) colore=1 else colore="gray"
  normpoints(q, type="l", lty=cont, col=colore)
  cont <- cont+1
  lableg <- c(lableg, "LP3")
  colleg <- c(colleg, colore)
 }
 legend("bottomright", legend=lableg, lty=c(1:cont), col=colleg, bty="n")
}


# ---------------------------------------------------- #

.lnML <- function(sample, dist="NORM") {
 # logarithm of the maximum likelihood
 # sample = sample
 # dist = distribution ("NORM","LN","EV1","EV2","GEV","GAM","LP3")
 n <- length(sample)
 if (any(c("LN","EV2","LP3")==dist)) sample <- log(sample)
 b <- sort(sample)
 if (any(is.nan(b))) {
  F <- NA
 }
 else {
  if (dist=="NORM") {
   T <- ML_estimation(b, dist="NORM")
   lnML <- sum(log(dnorm(b, T[1], T[2])))
  }
  else if (dist=="LN") {
   T <- ML_estimation(b, dist="NORM")
   lnML <- sum(log(dnorm(b, T[1], T[2])/exp(b)))
  }
  else if ((dist=="EV1")||(dist=="GUMBEL")) {
   T <- ML_estimation(b, dist="EV1")
   #lnML <- -.logLgumb(T, b)
   lnML <- sum(log(f.gumb(b, T[1], T[2])))
  }
  else if (dist=="EV2") {
   T <- ML_estimation(b, dist="EV1")
   lnML <- sum(log(f.gumb(b, T[1], T[2])/exp(b)))
  }
  else if (dist=="GEV") {
   T <- ML_estimation(b, dist="GEV")
   #lnML <- -.logLgev(T, b)
   if (abs(T[1] + T[2]/T[3] - b[n]) < 0.0001) {
    b <- b[-n]   # delete the maximum
    n2 <- n-1
   }
   else {
    n2 <- n
   }
   lnML <- n/n2 * sum(log(f.GEV(b, T[1], T[2], T[3])))
  }
  else if ((dist=="GAM")||(dist=="P3")) {
   T <- ML_estimation(b, dist="GAM")
   #lnML <- -.logLgam(T[1], b)
   if (abs(T[1] - b[1]) < 0.000001) {
    b <- b[-1]   # delete the minimum
    n2 <- n-1
   }
   else if (abs(T[1] - b[n]) < 0.000001) {
    b <- b[-n]   # delete the maximum
    n2 <- n-1
   }
   else {
    n2 <- n
   }
   lnML <- n/n2 * sum(log(f.gamma(b, T[1], T[2], T[3])))
  }
  else if (dist=="LP3") {
   T <- ML_estimation(b, dist="GAM")
   if (abs(T[1] - b[1]) < 0.000001) {
    b <- b[-1]   # delete the minimum
    n2 <- n-1
   }
   else if (abs(T[1] - b[n]) < 0.000001) {
    b <- b[-n]   # delete the maximum
    n2 <- n-1
   }
   else {
    n2 <- n
   }
   lnML <- n/n2 * sum(log(f.gamma(b, T[1], T[2], T[3])/exp(b)))
  }
  else stop(".lnML (sample, dist): distribution unknown")
 }
 return(lnML)
}


# ---------------------------------------------------- #

.AIC <- function(sample, dist="NORM") {
 # sample = sample
 # dist = distribution ("NORM","LN","EV1","EV2","GEV","GAM","LP3")
 if (any(c("NORM","LN","EV1","GUMBEL","EV2") == dist)) {
  p <- 2
 }
 else if (any(c("GEV","GAM","P3","LP3") == dist)) {
  p <- 3
 }
 else stop("AIC (sample, dist): distribution unknown")
 logL <- .lnML(sample, dist)
 AIC <- -2*logL + 2*p
 return(AIC)
}


# ---------------------------------------------------- #

.AICc <- function(sample, dist="NORM") {
 # sample = sample
 # dist = distribution ("NORM","LN","EV1","EV2","GEV","GAM","LP3")
 n <- length(sample)
 if (any(c("NORM","LN","EV1","GUMBEL","EV2") == dist)) {
  p <- 2
 }
 else if (any(c("GEV","GAM","P3","LP3") == dist)) {
  p <- 3
 }
 else stop("AIC (sample, dist): distribution unknown")
 logL <- .lnML(sample, dist)
 AICc <- -2*logL + 2*p*(n/(n-p-1))
 return(AICc)
}


# ---------------------------------------------------- #

.BIC <- function(sample, dist="NORM") {
 # sample = sample
 # dist = distribution ("NORM","LN","EV1","EV2","GEV","GAM","LP3")
 n <- length(sample)
 if (any(c("NORM","LN","EV1","GUMBEL","EV2") == dist)) {
  p <- 2
 }
 else if (any(c("GEV","GAM","P3","LP3") == dist)) {
  p <- 3
 }
 else stop("BIC (sample, dist): distribution unknown")
 logL <- .lnML(sample, dist)
 BIC <- -2*logL + log(n)*p
 return(BIC)
}


# ---------------------------------------------------- #

.ADC <- function(sample, dist="NORM") {
 # sample = sample
 # dist = distribution ("NORM","LN","EV1","EV2","GEV","GAM","LP3")
 if (any(c("LN","EV2","LP3")==dist)) sample <- log(sample)
 b <- sort(sample)
 n <- length(sample)
 eta0=0.851; beta0=0.116; csi0=0.0403; eps1=1.2; eps2=0.2
 if ((dist=="NORM")||(dist=="LN")) {
  T <- ML_estimation(b,dist="NORM")
  F <- .Fx(b,T,dist="NORM")
  eta1=1.147; beta1=0.229; csi1=0.167
  eta1corr <- eta1*(1+0.5/n); beta1corr <- beta1*(1-0.2/n); csi1corr <- csi1*(1+0.3/n)
 }
 else if ((dist=="EV1")||(dist=="GUMBEL")||(dist=="EV2")) {
  T <- ML_estimation(b,dist="EV1")
  F <- .Fx(b,T,dist="EV1")
  eta1=1.141; beta1=0.229; csi1=0.169
  eta1corr <- eta1*(1+0.5/n); beta1corr <- beta1*(1-0.2/n); csi1corr <- csi1*(1+0.1/n)
 }
 else if (dist=="GEV") {
  T <- ML_estimation(b,dist="GEV")
  F <- .Fx(b,T,dist="GEV")
  if (T[3]>0.5) T[3]=0.5
  eta1 <- 1.186*(1-0.04*T[3]-0.04*T[3]^2-0.01*T[3]^3)
  beta1 <- 0.189*(1+0.20*T[3]+0.37*T[3]^2+0.17*T[3]^3)
  csi1 <- 0.147*(1+0.13*T[3]+0.21*T[3]^2+0.09*T[3]^3)
  eta1corr <- eta1*(1-0.7/n+0.2/sqrt(n))
  beta1corr <- beta1*(1-1.8/n)
  csi1corr <- csi1*(1+0.9/n-0.2/sqrt(n))
 }
 else if ((dist=="GAM")||(dist=="P3")||(dist=="LP3")) {
  T <- ML_estimation(b,dist="GAM")
  if (abs(T[1] - b[1]) < 0.000001) {
   b <- b[-1]   # delete the minimum
   n <- n-1
  }
  else if (abs(T[1] - b[n]) < 0.000001) {
   b <- b[-n]   # delete the maximum
   n <- n-1
  }
  F <- .Fx(b,T,dist="GAM")
  F[F>0.99999999]=0.99999999
  F[F<0.00000001]=0.00000001
  if (T[3]<2) T[3]=2
  eta1 <- 1.194*(1-0.04*T[3]^(-1)-0.12*T[3]^(-2))
  beta1 <- 0.186*(1+0.34*T[3]^(-1)+0.30*T[3]^(-2))
  csi1 <- 0.145*(1+0.17*T[3]^(-1)+0.33*T[3]^(-2))
  eta1corr <- eta1*(1-1.8/n+0.1/sqrt(n)+0.5/(sqrt(n)*T[3]))
  beta1corr <- beta1*(1-0.5/n-0.3/sqrt(n)+0.3/(sqrt(n)*T[3]))
  csi1corr <- csi1*(1+2.0/n-0.3/sqrt(n)-0.4/(sqrt(n)*T[3]))
 }
 else stop("A2_GOFlaio(sample,T,dist,case): distribution unknown")
 A <- A2(F)
 if (A <= eps1*csi1corr) {
  ADC <- max((csi0 +beta0 *((eps1-1)*csi1corr/beta1corr)^(eta1corr/eta0))/((eps1-eps2)*csi1corr)*(A-eps2*csi1corr), 0.00001)
 }
 else {
  ADC <- csi0 + beta0*((A-csi1corr)/beta1corr)^(eta1corr/eta0)
 }
 if(any(is.nan(sample))) ADC <- NA
 return(ADC)
}

