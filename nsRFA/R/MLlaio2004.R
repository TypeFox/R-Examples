.Fx <- function (x, T, dist="NORM") {
 # x = sample
 # T = parameters (position, scale, shape, ...)
 # dist = distribution ("NORM","EV1","GEV","GAM")
 if (dist=="NORM") {
  pnorm(x,mean=T[1],sd=T[2])
 }
 else if ((dist=="EV1")||(dist=="GUMBEL")) {
  exp(-exp(-(x-T[1])/T[2]))
 }
 else if (dist=="GEV") {
  exp(-(1-((T[3]*(x-T[1]))/T[2]))^(1/T[3]))
 }
 else if ((dist=="GAM")||(dist=="P3")) {
  if (T[2]>0) {
   pgamma((x-T[1]), shape=T[3], scale=T[2])
   #pgamma((x-T[1])/T[2], shape=T[3])
  }
  else {
   1 - pgamma((T[1]-x), shape=T[3], scale=-T[2])
  }
 }
 else stop(".Fx(x, T, dist): distribution unknown")
}

# ------------------------------------------------------------------- #

moment_estimation <- function (x, dist="NORM") {
 # estimates the parameters of the distribution dist with the method of moments 
 # x = sample
 # dist = distribution ("NORM","EV1","GEV","GAM") 
 n <- length(x)
 mux <- mean(x)
 sdx <- sd(x)
 skx <- sum((x - mux)^3)/(n * sdx^3)
 if (dist=="NORM") {
  c(mux, sdx)
 }
 else if ((dist=="EV1")||(dist=="GUMBEL")) {
  muy <- 0.577216
  sdy <- pi/6^0.5
  c(mux-sdx*muy/sdy, sdx/sdy)
 }
 else if (dist=="GEV") {
  data(functionsLaio, envir=environment())
  #df.kgev <- read.table("kgev.dat") # df.kgev is a matrix with the skewness coefficient (first column) 
                                 # and the corresponding shape parameter of the GEV (second column)
  df.kgev  <- get("df.kgev", envir=environment())
  T3 <- approx(df.kgev[,1], df.kgev[,2], skx)$y
  T2 <- abs(T3)*sdx/(gamma(1+2*T3)-(gamma(1+T3))^2)^0.5
  T1 <- mux-T2/T3*(1-gamma(1+T3))
  c(T1,T2,T3)
 }
 else if ((dist=="GAM")||(dist=="P3")) {
  T3 <- 4/(skx^2)
  T2 <- sign(T3)*sdx/T3^0.5
  T1 <- mux-T2*T3
  c(T1,T2,T3)
 }
 else stop("moment_estimation(x, dist): distribution unknown")
}


# ------------------------------------------------------------------- #

.logLgumb <- function (T, x) {
 # .logLgumb is the negative log-likelihood function for the ev1 distribution
 # (T is the parameter vector and x is the sample)
 -sum(-(x-T[1])/T[2]-exp(-(x-T[1])/T[2])-log(T[2]))
}


# ------------------------------------------------------------------- #

.logLgev <- function (T, x) {
 # .logLgev is the negative log-likelihood function for the GEV distribution,
 # (T is the vector of parameters and x is the sample)
 n <- length(x)
 if (T[3]>1) warning(".logLgev(T,x): maximum likelihood estimator may not exist for shape parameter >1")
 if ((-0.0000001<T[3])&&(T[3]<0.0000001)) { # gumbel distribution for T[3]=0
  y <- (x-T[1])/T[2]
 }
 else {
  y <- -1/T[3]*log(1-T[3]*(x-T[1])/T[2])
 }
 n*log(T[2])+(1-T[3])*sum(y)+sum(exp(-y))
}


# ------------------------------------------------------------------- #

.logLgam <- function (T1,x) {
 # .logLgam is the negative log-likelihood function for the gamma distribution,
 # (T1 is the position parameter and x is the sample)
 n <- length(x)
 T2 <- 1/n*sum(x-T1)-n*(sum((x-T1)^(-1)))^(-1)  # scale parameter, Johnson et al. [1994], eq. 17.45
 T3 <- 1/n/T2*sum(x-T1)  # shape parameter, Johnson et al. [1994], eq. 17.45
 n*log(abs(T2))+n*lgamma(T3)-(T3-1)*sum(log((x-T1)/T2))+1/T2*sum(x-T1)
}


# ------------------------------------------------------------------- #

ML_estimation <- function (x, dist="NORM") {
 # estimates the parameters of the distribution dist with the method of maximum likelihood
 # x = sample
 # dist = distribution ("NORM","EV1","GEV","GAM")
 n <- length(x)
 Tm <- moment_estimation(x,dist)
 if (dist=="NORM") {
  T1 <- Tm[1]
  T2 <- Tm[2]*((n-1)/n)^0.5 # because ML estimotors require the variance 
                            # to be normilized by n (not n-1) 
  c(T1,T2)
 }
 else if ((dist=="EV1")||(dist=="GUMBEL")) {
  suppressWarnings(optim(par=Tm, fn=.logLgumb, x=x)$par)
 }
 else if (dist=="GEV") {
  suppressWarnings(optim(par=Tm, fn=.logLgev, x=x)$par)
 }
 else if ((dist=="GAM")||(dist=="P3")) {
  skx <- sum((x - mean(x))^3)/(length(x) * sd(x)^3)
  if (skx>0) { # this is for detecting if the distribution is lower or upper bounded
   Tm1 <- min(x)-1; a1 <- 1; a2 <- min(x)*(1-sign(min(x))*10^(-4))
  }
  else if (skx<0) {
   Tm1 <- max(x)+1; a1 <- -1; a2 <- -max(x)*(1+sign(max(x))*10^(-4))
  }
  T1 <- suppressWarnings(optim(par=Tm1, fn=.logLgam, x=x)$par)
  T2 <- 1/n*sum(x-T1)-n*(sum((x-T1)^(-1)))^(-1)    # scale parameter, Johnson et al. [1994], eq. 17.45
  T3 <- 1/n/T2*sum(x-T1)   # shape parameter, Johnson et al. [1994], eq. 17.45
  T4 <- 0
  if ((T3<2)&&(T2>0)) {   # non-regular estimation, see Smith [1985]
   data(functionsLaio, envir=environment())
   #df.polig <- read.table("polig.dat")  # Funzione poligamma
   df.polig <- get("df.polig", envir=environment())
   T4 <- 1
   b <- sort(x)
   T1 <- b[1]    # position parameter= first order statistic
   z <- b[2:n]-T1
   T3 <- approx(df.polig[,2], df.polig[,1], sum(log(z))/(n-1)-log(mean(z)))$y  # shape parameter, Johnson et al. [1994], eq. 17.48
   T2 <- mean(z)/T3
  }
  else if ((T3<2)&&(T2<0)) {
   data(functionsLaio, envir=environment())
   #df.polig <- read.table("polig.dat")
   df.polig <- get("df.polig", envir=environment())
   T4 <- 1
   b <- sort(x)
   T1 <- b[n]    # position parameter= n-th order statistic
   z <- b[1:(n-1)]-T1
   T3 <- approx(df.polig[,2], df.polig[,1], sum(log(abs(z)))/(n-1)-log(mean(abs(z))))$y  # shape parameter, Johnson et al. [1994], eq. 17.45
   T2 <- mean(z)/T3
  }
  c(T1,T2,T3)
 }
 else stop("ML_estimation(x, dist): distribution unknown")
}


# ----------------------------------------------------------------------------------- #

.sample_generator <- function (n, T, dist="NORM") {
 # n = sample length
 # T = parameters (position, scale, shape, ...)
 # dist = distribution ("NORM","EV1","GEV","GAM","EXP")
 if (dist=="NORM") {
  rnorm(n,mean=T[1],sd=T[2])
 }
 else if ((dist=="EV1")||(dist=="GUMBEL")) {
  q <- runif(n, min=0.0000000001, max=0.9999999999)
  T[1] - T[2]*log(log(1/q))
 }
 else if (dist=="GEV") {
  q <- runif(n, min=0.0000000001, max=0.9999999999)
  (1-(log(1/q))^T[3])*T[2]/T[3]+T[1]
 }
 else if ((dist=="GAM")||(dist=="P3")) {
  if (T[2]>0) {
   T[1] + rgamma(n, shape=T[3], scale=T[2])
  }
  else {
   T[1] - rgamma(n, shape=T[3], scale=-T[2])
  }
 }
 else if (dist=="EXP") {
  q <- runif(n, min=0.0000000001, max=0.9999999999)
  T[1] - T[2]*log(1-q)
 }
 else stop(".sample_generator(n, T, dist): distribution unknown")
}


