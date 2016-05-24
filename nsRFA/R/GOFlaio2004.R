A2 <- function (F) {
 # F = cumulative distribution function
 if (any(F<=0)||any(F>=1)) stop("A2(F): must be 0<F<1")
 n <- length(F)
 -n-(1/n)*sum(seq(1,2*n-1,by=2)*log(F) + seq(2*n-1,1,by=-2)*log(1-F))
}

# ------------------------------------------------------------------ #

W2 <- function (F) {
 # F = cumulative distribution function
 if (any(F<=0)||any(F>=1)) stop("A2(F): must be 0<F<1")
 n <- length(F)
 sum((F - seq(1,2*n-1,by=2)/(2*n))^2) + 1/(12*n)
}


# ----------------------------------------------------------------------------------------- #

fw2 <- function (w) {
 # approximation of the probability distribution of w (first 2 terms) when H0 is true (Anderson-Darling, 1952)
 if (w < 1.2) {
  ((exp(-(1/16)/w)*besselK((1/16)/w,1/4)+1.11803*exp(-(25/16)/w)*besselK((25/16)/w,1/4))/(w^0.5))/pi
 }
 else {
  1-10^(-2.2*w - 0.4)
 }
}


# ----------------------------------------------------------------------------------------- #

A2_GOFlaio <- function (x, dist="NORM") {
 # x = sample
 # dist = distribution ("NORM","LN","EV1","EV2","GEV","GAM","LP3")
 if (any(c("LN","EV2","LP3")==dist)) x <- log(x)
 b <- sort(x)
 n <- length(x)
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
 else stop("A2_GOFlaio(x,T,dist,case): distribution unknown")
 A <- A2(F)
 if (A <= eps1*csi1corr) {
  w <- max((csi0+beta0*((eps1-1)*csi1corr/beta1corr)^(eta1corr/eta0))/((eps1-eps2)*csi1corr)*(A-eps2*csi1corr),0.00001)
 }
 else {
  w <- csi0+beta0*((A-csi1corr)/beta1corr)^(eta1corr/eta0)
 }
 pA <- fw2(w)
 output <- c(A,pA)
 names(output) <- c("A2","p(A2)")
 return(output)
}


# ------------------------------------------------------------------------------------ #

.typeIerrorA2_GOFlaio <- function (n, T, alfa=0.05, dist="NORM", Nsim=1000) {
 # n = samples length
 # T = parameters (position, scale, shape, ...) only for case=0
 # alfa = significance level of the test
 # dist = distribution ("NORM","LN","EV1","EV2","GEV","GAM","LP3")
 # Nsim = number of repliques
 A <- rep(NA,Nsim)
 if (dist=="NORM") {
  for (i in 1:Nsim) {
   x <- .sample_generator(n,T,dist="NORM")
   A[i] <- A2_GOFlaio(x,dist="NORM")[2] 
  }
 }
 else if (dist=="LN") {
  for (i in 1:Nsim) {
   x <- .sample_generator(n,T,dist="NORM")
   x <- log(x)
   A[i] <- A2_GOFlaio(x,dist="LN")[2]
  }
 }
 else if ((dist=="EV1")||(dist=="GUMBEL")) {
  for (i in 1:Nsim) {
   x <- .sample_generator(n,T,dist="EV1")
   A[i] <- A2_GOFlaio(x,dist="EV1")[2]
  }
 }
 else if (dist=="EV2") {
  for (i in 1:Nsim) {
   x <- .sample_generator(n,T,dist="EV1")
   x <- log(x)
   A[i] <- A2_GOFlaio(x,dist="EV2")[2]
  }
 }
 else if (dist=="GEV") {
  for (i in 1:Nsim) {
   x <- .sample_generator(n,T,dist="GEV")
   A[i] <- A2_GOFlaio(x,dist="GEV")[2]
  }
 }
 else if ((dist=="GAM")||(dist=="P3")) {
  for (i in 1:Nsim) {
   x <- .sample_generator(n,T,dist="GAM")
   A[i] <- A2_GOFlaio(x,dist="GAM")[2]
  }
 }
 else if (dist=="LP3") {
  for (i in 1:Nsim) {
   x <- .sample_generator(n,T,dist="GAM")
   x <- log(x)
   A[i] <- A2_GOFlaio(x,dist="LP3")[2]
  }
 }
 else stop("typeIerror(n,T,dist,case,Nsim): distribution unknown")
 1 - ecdf(A)(1-alfa)
}

 # examples:
 #  > .typeIerrorA2_GOFlaio(100,c(10,3),0.1,dist="NORM",Nsim=10000)
 #  [1] 0.0979
 #  > .typeIerrorA2_GOFlaio(100,c(10,3),0.1,dist="EV1",Nsim=10000)
 #  [1] 0.0921
 #  > .typeIerrorA2_GOFlaio(100,c(10,3,0.1),0.1,dist="GEV",Nsim=1000)
 #  [1] 0.092
 # spesso con la GEV mi da l'errore:
 #  Error in optim(par = Tm, fn = .logLgev, x = x) :
 #        function cannot be evaluated at initial parameters
 #  > .typeIerrorA2_GOFlaio(100,c(-10,0.3,70),0.1,dist="GAM",Nsim=10000)
 #  [1] 0.0913
 #  > .typeIerrorA2_GOFlaio(100,c(20,1,70),0.1,dist="GAM",Nsim=10000)
 #  [1] 0.0949
 #  > .typeIerrorA2_GOFlaio(100,c(0.2800572,0.1237543,5.817518),0.1,dist="GAM",Nsim=10000)
 #  [1] 0.0863
 #  > .typeIerrorA2_GOFlaio(100,c(226.714,130.4902,2.216652),0.1,dist="GAM",Nsim=10000)
 #  [1] 0.1993 ????   0.1951
 #  > .typeIerrorA2_GOFlaio(100,c(167.157813,99.583519,3.502665),0.1,dist="GAM",Nsim=10000)
 #  [1] 0.1026
 #  > .typeIerrorA2_GOFlaio(100,c(0.09,0.08,10),0.1,dist="GAM",Nsim=10000)
 #  [1] 0.0953
 #  > .typeIerrorA2_GOFlaio(100,c(100,3),0.1,dist="LN",Nsim=10000)
 #  [1] 0.1120
 #  > .typeIerrorA2_GOFlaio(100,c(100,3),0.1,dist="EV2",Nsim=10000)
 #  [1] 0.1149
 #  > .typeIerrorA2_GOFlaio(100,c(-10,1,70),0.1,dist="LP3",Nsim=10000)
 #  [1] 0.1043
 #  > .typeIerrorA2_GOFlaio(100,c(100,3),0.1,dist="EV1",Nsim=10000)
 #  [1] 0.1003

