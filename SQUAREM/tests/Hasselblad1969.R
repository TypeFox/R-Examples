require("SQUAREM")  # master funcion and the 4 different SQUAREM algorithms
require("setRNG")   # needed to reproduce the same results

# data for Poisson mixture estimation  from Hasselblad (JASA 1969)
poissmix.dat <- data.frame(death=0:9, freq=c(162,267,271,185,111,61,27,8,3,1))

y <- poissmix.dat$freq
tol <- 1.e-08

# generate a random initial guess for 3 parameters
setRNG(list(kind="Wichmann-Hill", normal.kind="Box-Muller", seed=123))
p0 <- c(runif(1),runif(2,0,4))    

# fixptfn = function that computes a single update of fixed-point iteration
# objfn = underlying objective function that needs to be maximized
#

# EM algorithm

poissmix.em <- function(p,y) {
# The fixed point mapping giving a single E and M step of the EM algorithm
# 
pnew <- rep(NA,3)
i <- 0:(length(y)-1)
zi <- p[1]*exp(-p[2])*p[2]^i / (p[1]*exp(-p[2])*p[2]^i + (1 - p[1])*exp(-p[3])*p[3]^i)
pnew[1] <- sum(y*zi)/sum(y)
pnew[2] <- sum(y*i*zi)/sum(y*zi)
pnew[3] <- sum(y*i*(1-zi))/sum(y*(1-zi))
p <- pnew
return(pnew)
}

poissmix.loglik <- function(p,y) {
# Objective function whose local minimum is a fixed point \
# negative log-likelihood of binary poisson mixture
i <- 0:(length(y)-1)
loglik <- y*log(p[1]*exp(-p[2])*p[2]^i/exp(lgamma(i+1)) + 
		(1 - p[1])*exp(-p[3])*p[3]^i/exp(lgamma(i+1)))
return ( -sum(loglik) )
}

good.tol <- 1.e-05

good <-  c(0.6401146029910, 2.6634043566619, 1.2560951012662)

  pf1 <- fpiter(p=p0, y=y, fixptfn=poissmix.em, objfn=poissmix.loglik, control=list(tol=tol))
  print(pf1$par, digits=16)
#  good <- c(0.6401136228187346, 2.6634055536337877, 1.2560968051580201)
  if (good.tol < max(abs(good - pf1$par))) stop("error in EM algorithm.")


# First-order SQUAREM algorithm with SqS3 method
  pf2 <- squarem(par=p0, y=y, fixptfn=poissmix.em, objfn=poissmix.loglik,   
    control=list(tol=tol))
  print(pf2$par, digits=16)
#[1] 0.640114475404004 2.663404513995809 1.256095320343123 linux Ubuntu 10.04 Intel core 2 duo/Dell
#[1] 0.640114475330512 2.663404514085158 1.256095320471585 linux Ubuntu 10.04 Intel Centrino/Toshiba

#  good <- c(0.640114475404004, 2.663404513995809, 1.256095320343123)
  if (good.tol < max(abs(good - pf2$par))) stop("error in First-order SQUAREM algorithm with SqS3 method.")

# First-order SQUAREM algorithm with SqS2 method
  pf3 <- squarem(par=p0, y=y, fixptfn=poissmix.em, objfn=poissmix.loglik, 
    control=list(method=2, tol=tol))
  print(pf3$par, digits=16)
#[1] 0.640115228256831  2.663403593105286  1.256094014312148 linux Ubuntu 10.04 Intel core 2 duo/Dell
#[1] 0.6401152283764715 2.6634035929591846 1.2560940141041679linux Ubuntu 10.04 Intel Centrino/Toshiba

#  good <- c(0.640115228256831, 2.663403593105286, 1.256094014312148)
  if (good.tol < max(abs(good - pf3$par))) stop("error in First-order SQUAREM algorithm with SqS2 method.")

# First-order SQUAREM algorithm with SqS3 method; non-monotone 
# Note: the objective function is not evaluated when objfn.inc = Inf 
  pf4 <- squarem(par=p0,y=y, fixptfn=poissmix.em, 
    control=list(tol=tol, objfn.inc=Inf))
  print(pf4$par, digits=16)
#[1] 0.640114475404004 2.663404513995809 1.256095320343123 linux Ubuntu 10.04 Intel core 2 duo/Dell
#[1] 0.640114475330512 2.663404514085158 1.256095320471585 linux Ubuntu 10.04 Intel Centrino/Toshiba

#  good <- c(0.640114475404004, 2.663404513995809, 1.256095320343123)
  if (good.tol < max(abs(good - pf4$par))) stop("error in First-order SQUAREM algorithm with SqS3 method; non-monotone.")

# First-order SQUAREM algorithm with SqS3 method; objective function is not specified
  pf5 <- squarem(par=p0,y=y, fixptfn=poissmix.em, control=list(tol=tol, kr=0.1))
  print(pf5$par, digits=16)
#[1] 0.640114475404004 2.663404513995809 1.256095320343123 linux Ubuntu 10.04 Intel core 2 duo/Dell
#[1] 0.640114475330512 2.663404514085158 1.256095320471585 linux Ubuntu 10.04 Intel Centrino/Toshiba

#  good <- c(0.640114475404004, 2.663404513995809, 1.256095320343123)
  if (good.tol < max(abs(good - pf5$par))) stop("error in First-order SQUAREM algorithm with SqS3 (and no objective function).")

# Second-order (K=2) SQUAREM algorithm with SqRRE 
  pf6 <- squarem(par=p0, y=y, fixptfn=poissmix.em, objfn=poissmix.loglik,
    control=list (K=2, tol=tol))
  print(pf6$par, digits=16)
#[1] 0.640114603181059 2.663404356560964 1.256095100702419 linux Ubuntu 10.04 Intel core 2 duo/Dell
#[1] 0.640114603175754 2.663404356262396 1.256095101254215 linux Ubuntu 10.04 Intel Centrino/Toshiba

#  good <- c(0.640114603181059, 2.663404356560964, 1.256095100702419)
  if (good.tol < max(abs(good - pf6$par))) stop("error in Second-order (K=2) SQUAREM algorithm with SqRRE.")

# Second-order SQUAREM algorithm with SqRRE; objective function is not specified
  pf7 <- squarem(par=p0, y=y, fixptfn=poissmix.em, control=list(K=2, tol=tol))
  print(pf7$par, digits=16)
#[1] 0.640114603108274  2.663404357034439  1.256095100144889  linux Ubuntu 10.04 Intel core 2 duo/Dell
#[1] 0.6401146029029904 2.6634043571424155 1.2560951007555829 linux Ubuntu 10.04 Intel Centrino/Toshiba

#  good <- c(0.640114603108274, 2.663404357034439, 1.256095100144889)
  if (good.tol < max(abs(good - pf7$par))) stop("error in Second-order SQUAREM algorithm with SqRRE (and no objective function).")

# Number of fixed-point evaluations
  c(pf1$fpeval, pf2$fpeval, pf3$fpeval, pf4$fpeval, pf5$fpeval, pf6$fpeval, pf7$fpeval)

# Number of objective function evaluations
  c(pf1$objfeval, pf2$objfeval, pf3$objfeval, pf4$objfeval, pf5$objfeval, pf6$objfeval, pf7$objfeval)

# Comparison of converged parameter estimates
  par.mat <- rbind(pf1$par, pf2$par, pf3$par, pf4$par, pf5$par, pf6$par, pf7$par)
  par.mat

