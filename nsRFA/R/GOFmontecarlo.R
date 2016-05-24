 
# gofP3test <- function (x,Nsim=1000) {
# 
#  # Monte-Carlo procedure
#  # INPUT:
#  # x = sample
#  # Nsim = number of generations
#  x <- sort(x)
#  n <- length(x)
#  Lmom.x <- Lmoments(x)
# 
#  lambda1o <- Lmom.x["l1"]
#  lambda2o <- Lmom.x["l2"]
#  tau3o <- Lmom.x["lca"]
#  if ((abs(tau3o) > 0)&&(abs(tau3o) < 1/3)) {
#   z <- 3*pi*tau3o^2
#   alfa0 <- (1 + 0.2906*z)/(z + 0.1882*z^2 + 0.0442*z^3)
#  }
#  else if ((abs(tau3o) >= 1/3)&&(abs(tau3o) < 1)) {
#   z <- 1 - abs(tau3o)
#   alfa0 <- (0.36067*z - 0.59567*z^2 + 0.25361*z^3)/(1 - 2.78861*z + 2.56096*z^2 - 0.77045*z^3)
#  }
# 
#  if (alfa0<100) {
#   sigma <- lambda2o*pi^(0.5) * alfa0^(0.5) * gamma(alfa0)/gamma(alfa0 + 0.5)
#   beta0 <- 0.5*sigma*abs(2*alfa0^(-0.5))
#   if(beta0>0) {
#    xi0 <- lambda1o - alfa0*beta0
#    #F <- pgamma((x - xi0)/beta0, alfa0)
#    F <- pgamma(x-xi0,shape=alfa0,scale=beta0)
#   }
#   else {
#    xi0 <- lambda1o + alfa0*beta0
#    #F <- 1 - pgamma((xi0 - x)/beta0, alfa0)
#    F <- 1 - pgamma(xi0-x,shape=alfa0,scale=-beta0)
#   }
#  }
#  else {
#   sigma0 <- sqrt(pi)*lambda2o/(1-1/(8*alfa0)+1/(128*alfa0^2))
#   mu0 <- lambda1o
#   beta0 <- NA
#   F <- pnorm(x,mu0,sigma0)
#  }
#  F[F<=0] <- 0.00000001
#  F[F>=1] <- 0.99999999
#  A2 <- -n-(1/n)*sum((seq(1,2*n-1,by=2))*log(F) + (seq(2*n-1,1,by=-2))*log(1-F))
# 
#  A2s <- rep(NA,Nsim)
#  for (i in 1:Nsim) {
#   #F.sim <- sort(runif(n, min=0.0000000001, max=0.9999999999))
#   #if (is.na(beta0)) {
#   # x.sim <- qnorm(F.sim,mu0,sigma0)
#   #}
#   #else if (beta0 >= 0) {
#   #  x.st <- qgamma(F.sim, alfa0)
#   #  x.sim <- x.st*beta0 + xi0
#   #}
#   #else if (beta0 < 0) {
#   #  x.st <- qgamma(F.sim, alfa0)
#   #  x.sim <- xi0 - x.st*beta0
#   #}
#   if (is.na(beta0)) {
#    x.sim <- rnorm(n,mu0,sigma0)
#   }
#   else if (beta0 > 0) {
#    x.sim <- xi0 + rgamma(n,shape=alfa0,scale=beta0)
#   }
#   else {
#    x.sim <- xi0 - rgamma(n,shape=alfa0,scale=-beta0)
#   }
#   Lmom.x.sim <- Lmoments(x.sim)
#   lambda1 <- Lmom.x.sim["l1"]
#   lambda2 <- Lmom.x.sim["l2"]
#   tau3 <- Lmom.x.sim["lca"]
#   #cat("\n\nF.sim\n")
#   #print(F.sim)
#   #cat("\nL-moments\n")
#   #print(c(lambda1,lambda2,tau3))
#   #cat("\nx.st\n")
#   #print(x.st)
#   #cat("\nxi0, alfa0, beta0\n")
#   #print(c(xi0,alfa0,beta0))
#   if ((abs(tau3) > 0)&&(abs(tau3) < 1/3)) {
#    z <- 3*pi*tau3^2
#    alfa <- (1 + 0.2906*z)/(z + 0.1882*z^2 + 0.0442*z^3)
#   }
#   else if ((abs(tau3) >= 1/3)&&(abs(tau3) < 1)) {
#    z <- 1 - abs(tau3)
#    alfa <- (0.36067*z - 0.59567*z^2 + 0.25361*z^3)/(1 - 2.78861*z + 2.56096*z^2 - 0.77045*z^3)
#   }
# 
#   if (alfa<100) {
#    sigma <- lambda2*pi^(0.5) * alfa^(0.5) * gamma(alfa)/gamma(alfa + 0.5)
#    beta <- 0.5*sigma*abs(2*alfa^(-0.5))
#    if(beta>0) {
#     xi <- lambda1 - alfa*beta
#     #F <- pgamma((x - xi)/beta, alfa)
#     F <- pgamma(x-xi,shape=alfa,scale=beta)
#    }
#    else {
#     xi <- lambda1 + alfa*beta
#     #F <- 1 - pgamma((xi - x)/beta, alfa)
#     F <- 1 - pgamma(xi-x,shape=alfa,scale=-beta)
#    }
#   }
#   else {
#    sigma <- sqrt(pi)*lambda2/(1-1/(8*alfa)+1/(128*alfa^2))
#    F <- pnorm(x,lambda1,sigma)
#   }
#   F[F<0.00000001] <- 0.00000001
#   F[F>0.99999999] <- 0.99999999
#   A2s[i] <- -n-(1/n)*sum((seq(1,2*n-1,by=2))*log(F) + (seq(2*n-1,1,by=-2))*log(1-F))
#  }
# 
#  ecdfA2s <- ecdf(A2s)
#  probabilita <- ecdfA2s(A2)
#  output <- signif(c(A2, probabilita),4)
#  names(output) <- c("A2","P")
# 
#  return(output)
# }
 
 
# -------------------------------------------------------------------------------------- #
 
gofP3test <- function (x,Nsim=1000) {

 # Monte-Carlo procedure
 # INPUT:
 # x = sample
 # Nsim = number of generations
 x <- sort(x)
 n <- length(x)
 Lmom.x <- Lmoments(x)

 par <- par.gamma(Lmom.x["l1"],Lmom.x["l2"],Lmom.x["lca"])
 if(is.na(par$beta)) F <- pnorm(x,par$mu,par$sigma)
 else F <- F.gamma(x,par$xi,par$beta,par$alfa)
 F[F<0.00000001] <- 0.00000001
 F[F>0.99999999] <- 0.99999999
 A2 <- -n-(1/n)*sum((seq(1,2*n-1,by=2))*log(F) + (seq(2*n-1,1,by=-2))*log(1-F))

 A2s <- rep(NA,Nsim)
 for (i in 1:Nsim) {
  if(!is.na(par$beta)) x.sim <- rand.gamma(n,par$xi,par$beta,par$alfa)
  else x.sim <- rnorm(n,par$mu,par$sigma)
  x.sim <- sort(x.sim)
  Lmom.xsim <- Lmoments(x.sim)
  par.sim <- par.gamma(Lmom.xsim["l1"],Lmom.xsim["l2"],Lmom.xsim["lca"])
  if(is.na(par.sim$beta)) F <- pnorm(x.sim,par.sim$mu,par.sim$sigma)
  else F <- F.gamma(x.sim,par.sim$xi,par.sim$beta,par.sim$alfa)
  F[F<0.00000001] <- 0.00000001
  F[F>0.99999999] <- 0.99999999
  A2s[i] <- -n-(1/n)*sum((seq(1,2*n-1,by=2))*log(F) + (seq(2*n-1,1,by=-2))*log(1-F))
 }

 ecdfA2s <- ecdf(A2s)
 probabilita <- ecdfA2s(A2)
 output <- signif(c(A2, probabilita),4)
 names(output) <- c("A2","P")

 return(output)
}


# ----------------------------------------------------------------------------- #

gofNORMtest <- function(x) {

 # Francesco Laio

 fw2 <- function (x) {
  if (x < 1.2) {
   fw2 <- ((exp(-(1/16)/x)*besselK((1/16)/x,1/4)+1.11803*exp(-(25/16)/x)*besselK((25/16)/x,1/4))/(x^0.5))/pi
  }
  else fw2 <- 1
  return(fw2)
 }

 b <- sort(x)
 n <-length(x)
 c0=0.851; b0=0.116; csi0=0.0403; eps1=1.2; eps2=0.2

 Mgaus <- c(mean(b),sum((b-mean(b))^2 / n)^0.5)
 Fgaus <- pmax(pmin(pnorm(b,Mgaus[1],Mgaus[2]),0.999999999),0.00000001)

 # Anderson-Darling

 c1=1.147; b1=0.229; csi1=0.167
 c1corr <- c1*(1+0.5/n); b1corr <- b1*(1-0.2/n); csi1corr <- csi1*(1+0.3/n)

 Agaus <- -n-sum((2*(1:n)-1)*log(Fgaus)+(2*n+1-2*(1:n))*log(1-Fgaus))/n

 if (Agaus <= eps1*csi1corr) {
   z3 <- max((csi0+b0*((eps1-1)*csi1corr/b1corr)^(c1corr/c0))/((eps1-eps2)*csi1corr)*(Agaus-eps2*csi1corr),0.00001)
 }
 else {
   z3 <- csi0+b0*((Agaus-csi1corr)/b1corr)^(c1corr/c0)
 }

 pA <- fw2(z3)
 
 output <- c(Agaus,pA); names(output) <- c("A2","P")
 return(output)
}


# -------------------------------------------------------------------------------------- #

gofGEVtest <- function (x,Nsim=1000) {

 # Monte-Carlo procedure
 # INPUT:
 # x = sample
 # Nsim = number of generations
 x <- sort(x)
 n <- length(x)
 Lmom.x <- Lmoments(x)

 par <- par.GEV(Lmom.x["l1"],Lmom.x["l2"],Lmom.x["lca"])
 F <- suppressWarnings(F.GEV(x,par$xi,par$alfa,par$k))
 F[F<0.00000001] <- 0.00000001
 F[F>0.99999999] <- 0.99999999
 A2 <- -n-(1/n)*sum((seq(1,2*n-1,by=2))*log(F) + (seq(2*n-1,1,by=-2))*log(1-F))

 A2s <- rep(NA,Nsim)
 for (i in 1:Nsim) {
  x.sim <- rand.GEV(n,par$xi,par$alfa,par$k)
  x.sim <- sort(x.sim)
  Lmom.xsim <- Lmoments(x.sim)
  par.sim <- par.GEV(Lmom.xsim["l1"],Lmom.xsim["l2"],Lmom.xsim["lca"])
  F <- suppressWarnings(F.GEV(x.sim,par.sim$xi,par.sim$alfa,par.sim$k))
  F[F<0.00000001] <- 0.00000001
  F[F>0.99999999] <- 0.99999999
  A2s[i] <- -n-(1/n)*sum((seq(1,2*n-1,by=2))*log(F) + (seq(2*n-1,1,by=-2))*log(1-F))
 }

 ecdfA2s <- ecdf(A2s)
 probabilita <- ecdfA2s(A2)
 output <- signif(c(A2, probabilita),4)
 names(output) <- c("A2","P")

 return(output)
}


# -------------------------------------------------------------------------------------- #

gofGENLOGIStest <- function (x,Nsim=1000) {

 # Monte-Carlo procedure
 # INPUT:
 # x = sample
 # Nsim = number of generations
 x <- sort(x)
 n <- length(x)
 Lmom.x <- Lmoments(x)

 par <- par.genlogis(Lmom.x["l1"],Lmom.x["l2"],Lmom.x["lca"])
 F <- suppressWarnings(F.genlogis(x,par$xi,par$alfa,par$k))
 F[F<0.00000001] <- 0.00000001
 F[F>0.99999999] <- 0.99999999
 A2 <- -n-(1/n)*sum((seq(1,2*n-1,by=2))*log(F) + (seq(2*n-1,1,by=-2))*log(1-F))

 A2s <- rep(NA,Nsim)
 for (i in 1:Nsim) {
  x.sim <- rand.genlogis(n,par$xi,par$alfa,par$k)
  x.sim <- sort(x.sim)
  Lmom.xsim <- Lmoments(x.sim)
  par.sim <- par.genlogis(Lmom.xsim["l1"],Lmom.xsim["l2"],Lmom.xsim["lca"])
  F <- suppressWarnings(F.genlogis(x.sim,par.sim$xi,par.sim$alfa,par.sim$k))
  F[F<0.00000001] <- 0.00000001
  F[F>0.99999999] <- 0.99999999
  A2s[i] <- -n-(1/n)*sum((seq(1,2*n-1,by=2))*log(F) + (seq(2*n-1,1,by=-2))*log(1-F))
 }

 ecdfA2s <- ecdf(A2s)
 probabilita <- ecdfA2s(A2)
 output <- signif(c(A2, probabilita),4)
 names(output) <- c("A2","P")

 return(output)
}


# -------------------------------------------------------------------------------------- #

gofGENPARtest <- function (x,Nsim=1000) {

 # Monte-Carlo procedure
 # INPUT:
 # x = sample
 # Nsim = number of generations
 x <- sort(x)
 n <- length(x)
 Lmom.x <- Lmoments(x)

 par <- par.genpar(Lmom.x["l1"],Lmom.x["l2"],Lmom.x["lca"])
 F <- suppressWarnings(F.genpar(x,par$xi,par$alfa,par$k))
 F[F<0.00000001] <- 0.00000001
 F[F>0.99999999] <- 0.99999999
 F[is.nan(F)] <- 0.99999999
 A2 <- -n-(1/n)*sum((seq(1,2*n-1,by=2))*log(F) + (seq(2*n-1,1,by=-2))*log(1-F))

 A2s <- rep(NA,Nsim)
 for (i in 1:Nsim) {
  x.sim <- rand.genpar(n,par$xi,par$alfa,par$k)
  x.sim <- sort(x.sim)
  Lmom.xsim <- Lmoments(x.sim)
  par.sim <- par.genpar(Lmom.xsim["l1"],Lmom.xsim["l2"],Lmom.xsim["lca"])
  F <- suppressWarnings(F.genpar(x.sim,par.sim$xi,par.sim$alfa,par.sim$k))
  F[F<0.00000001] <- 0.00000001
  F[F>0.99999999] <- 0.99999999
  F[is.nan(F)] <- 0.99999999
  A2s[i] <- -n-(1/n)*sum((seq(1,2*n-1,by=2))*log(F) + (seq(2*n-1,1,by=-2))*log(1-F))
 }

 ecdfA2s <- ecdf(A2s)
 probabilita <- ecdfA2s(A2)
 output <- signif(c(A2, probabilita),4)
 names(output) <- c("A2","P")

 return(output)
}


# -------------------------------------------------------------------------------------- #

gofLOGNORMtest <- function (x,Nsim=1000) {

 # Monte-Carlo procedure
 # INPUT:
 # x = sample
 # Nsim = number of generations
 x <- sort(x)
 n <- length(x)
 Lmom.x <- Lmoments(x)

 par <- par.lognorm(Lmom.x["l1"],Lmom.x["l2"],Lmom.x["lca"])
 F <- suppressWarnings(F.lognorm(x,par$xi,par$alfa,par$k))
 F[F<0.00000001] <- 0.00000001
 F[F>0.99999999] <- 0.99999999
 A2 <- -n-(1/n)*sum((seq(1,2*n-1,by=2))*log(F) + (seq(2*n-1,1,by=-2))*log(1-F))

 A2s <- rep(NA,Nsim)
 for (i in 1:Nsim) {
  x.sim <- rand.lognorm(n,par$xi,par$alfa,par$k)
  x.sim <- sort(x.sim)
  Lmom.xsim <- Lmoments(x.sim)
  par.sim <- par.lognorm(Lmom.xsim["l1"],Lmom.xsim["l2"],Lmom.xsim["lca"])
  F <- suppressWarnings(F.lognorm(x.sim,par.sim$xi,par.sim$alfa,par.sim$k))
  F[F<0.00000001] <- 0.00000001
  F[F>0.99999999] <- 0.99999999
  A2s[i] <- -n-(1/n)*sum((seq(1,2*n-1,by=2))*log(F) + (seq(2*n-1,1,by=-2))*log(1-F))
 }

 ecdfA2s <- ecdf(A2s)
 probabilita <- ecdfA2s(A2)
 output <- signif(c(A2, probabilita),4)
 names(output) <- c("A2","P")

 return(output)
}


# -------------------------------------------------------------------------------------- #

gofEXPtest <- function (x,Nsim=1000) {

 # Monte-Carlo procedure
 # INPUT:
 # x = sample
 # Nsim = number of generations
 x <- sort(x)
 n <- length(x)
 Lmom.x <- Lmoments(x)

 par <- par.exp(Lmom.x["l1"],Lmom.x["l2"])
 F <- suppressWarnings(F.exp(x,par$xi,par$alfa))
 F[F<0.00000001] <- 0.00000001
 F[F>0.99999999] <- 0.99999999
 A2 <- -n-(1/n)*sum((seq(1,2*n-1,by=2))*log(F) + (seq(2*n-1,1,by=-2))*log(1-F))

 A2s <- rep(NA,Nsim)
 for (i in 1:Nsim) {
  x.sim <- rand.exp(n,par$xi,par$alfa)
  x.sim <- sort(x.sim)
  Lmom.xsim <- Lmoments(x.sim)
  par.sim <- par.exp(Lmom.xsim["l1"],Lmom.xsim["l2"])
  F <- suppressWarnings(F.exp(x.sim,par.sim$xi,par.sim$alfa))
  F[F<0.00000001] <- 0.00000001
  F[F>0.99999999] <- 0.99999999
  A2s[i] <- -n-(1/n)*sum((seq(1,2*n-1,by=2))*log(F) + (seq(2*n-1,1,by=-2))*log(1-F))
 }

 ecdfA2s <- ecdf(A2s)
 probabilita <- ecdfA2s(A2)
 output <- signif(c(A2, probabilita),4)
 names(output) <- c("A2","P")

 return(output)
}


# -------------------------------------------------------------------------------------- #

gofGUMBELtest <- function (x,Nsim=1000) {

 # Monte-Carlo procedure
 # INPUT:
 # x = sample
 # Nsim = number of generations
 x <- sort(x)
 n <- length(x)
 Lmom.x <- Lmoments(x)

 par <- par.gumb(Lmom.x["l1"],Lmom.x["l2"])
 F <- suppressWarnings(F.gumb(x,par$xi,par$alfa))
 F[F<0.00000001] <- 0.00000001
 F[F>0.99999999] <- 0.99999999
 A2 <- -n-(1/n)*sum((seq(1,2*n-1,by=2))*log(F) + (seq(2*n-1,1,by=-2))*log(1-F))

 A2s <- rep(NA,Nsim)
 for (i in 1:Nsim) {
  x.sim <- rand.gumb(n,par$xi,par$alfa)
  x.sim <- sort(x.sim)
  Lmom.xsim <- Lmoments(x.sim)
  par.sim <- par.gumb(Lmom.xsim["l1"],Lmom.xsim["l2"])
  F <- suppressWarnings(F.gumb(x.sim,par.sim$xi,par.sim$alfa))
  F[F<0.00000001] <- 0.00000001
  F[F>0.99999999] <- 0.99999999
  A2s[i] <- -n-(1/n)*sum((seq(1,2*n-1,by=2))*log(F) + (seq(2*n-1,1,by=-2))*log(1-F))
 }

 ecdfA2s <- ecdf(A2s)
 probabilita <- ecdfA2s(A2)
 output <- signif(c(A2, probabilita),4)
 names(output) <- c("A2","P")

 return(output)
}


# ---------------------------------------------------------------------------------- #

.test.GOFmontecarlo <- function (parameters, type="NORM", alfa=.05, n=30, N=100) {

  # testa se i test di goodness of fit funzionano correttamente
  # type = "NORM", "GENLOGIS", "GENPAR", "GEV", "LOGNORM", "P3", "EXP", "GUMBEL"
  # alfa = limite di significativita
  # n = lunghezza dei campioni
  # N = numero di ripetizioni

  ps <- rep(NA,N)
  if (type=="NORM") {
   for (i in 1:N) {
    x <- rnorm(n,parameters[1],parameters[2])
    ps[i] <- gofNORMtest(x)["P"]
   }
  }
  else if (type=="EXP") {
   for (i in 1:N) {
    x <- rand.exp(n,parameters[1],parameters[2])
    ps[i] <- gofEXPtest(x)["P"]
   }
  }
  else if (type=="GUMBEL") {
   for (i in 1:N) {
    x <- rand.gumb(n,parameters[1],parameters[2])
    ps[i] <- gofGUMBELtest(x)["P"]
   }
  }
  else if (type=="GENLOGIS") {
   for (i in 1:N) {
    x <- rand.genlogis(n,parameters[1],parameters[2],parameters[3])
    ps[i] <- gofGENLOGIStest(x)["P"]
   }
  }
  else if (type=="GENPAR") {
   for (i in 1:N) {
    x <- rand.genpar(n,parameters[1],parameters[2],parameters[3])
    ps[i] <- gofGENPARtest(x)["P"]
   }
  }
  else if (type=="GEV") { 
   for (i in 1:N) {
    x <- rand.GEV(n,parameters[1],parameters[2],parameters[3])
    ps[i] <- gofGEVtest(x)["P"]
   }
  }
  else if ((type=="LOGNORM")||(type=="LN3")) {
   for (i in 1:N) {
    x <- rand.lognorm(n,parameters[1],parameters[2],parameters[3])
    ps[i] <- gofLOGNORMtest(x)["P"]
   }
  }
  else if ((type=="P3")||(type=="GAM")) {
   for (i in 1:N) {
    x <- rand.gamma(n,parameters[1],parameters[2],parameters[3])
    ps[i] <- gofP3test(x)["P"]
   }
  }

  #typeIer <- sum(ps>(1-alfa))/N 
  typeIer <- 1 - ecdf(ps)(1-alfa)

  return(typeIer)
}

# .test.GOFmontecarlo(c(0,1),type="NORM",alfa=.1,N=10000) = 0.0976
# .test.GOFmontecarlo(c(10,1),type="NORM",alfa=.1,N=10000) = 0.0963
# .test.GOFmontecarlo(c(10,100),type="NORM",alfa=.1,N=10000) = 0.0953
# x <- annualflows[annualflows[1]==30,3]
# lmom <- Lmoments(x)
# param <- par.genlogis(lmom["l1"],lmom["l2"],lmom["lca"])
# .test.GOFmontecarlo(c(param$xi,param$alfa,param$k),type="GENLOGIS",n=100,alfa=.1,N=1000)  # = 0.137 
# param <- par.genpar(lmom["l1"],lmom["l2"],lmom["lca"])
# .test.GOFmontecarlo(c(param$xi,param$alfa,param$k),type="GENPAR",n=100,alfa=.1,N=1000)  # = 0.091 
# param <- par.GEV(lmom["l1"],lmom["l2"],lmom["lca"])
# .test.GOFmontecarlo(c(param$xi,param$alfa,param$k),type="GEV",n=100,alfa=.1,N=1000)  # = 0.099 
# param <- par.lognorm(lmom["l1"],lmom["l2"],lmom["lca"])
# .test.GOFmontecarlo(c(param$xi,param$alfa,param$k),type="LOGNORM",n=100,alfa=.1,N=1000)  # = 0.122
# param <- par.gamma(lmom["l1"],lmom["l2"],lmom["lca"])
# .test.GOFmontecarlo(c(param$xi,param$beta,param$alfa),type="P3",n=100,alfa=.1,N=1000)   # = 0.043  0.044
# .test.GOFmontecarlo(c(20,1,70),type="P3",n=100,alfa=.1,N=1000)   # = 0.103 
# .test.GOFmontecarlo(c(167.157813,99.583519,3.502665),type="P3",n=100,alfa=.1,N=1000)   # = 0.051
# .test.GOFmontecarlo(c(226.714,130.4902,2.216652),type="P3",n=100,alfa=.1,N=1000)   # = 0.032
# param <- par.exp(lmom["l1"],lmom["l2"])
# .test.GOFmontecarlo(c(param$xi,param$alfa),type="EXP",n=100,alfa=.1,N=1000)   # = 0.109
# param <- par.gumb(lmom["l1"],lmom["l2"])
# .test.GOFmontecarlo(c(param$xi,param$alfa),type="GUMBEL",n=100,alfa=.1,N=1000)   # = 0.077  0.113 

