######################################################################
# Compute actual coverage of an (1-\alpha)*100% confidence interval
# for a binomial proportion
#
# Params:
#  ci.fun - the binom.confint style function to compute the confidence
#           interval
#  n      - size of the binomial distribution
#  alpha  - significance level of the confidence interval
#  p.grid - vector of proportions where to evaluate. The function
#           automatically choses a grid of values where the
#           minimum coverage is to be found, but for plotting
#           purposes one might be interested in selecting more values for p
#           In case p.grid is not NULL the two set of values are merged.
# interval- In some situations one might be interesting in computing
#           coverage only for a restricted interval of c(0,1)
#  ...    - further arguments passed on to the ci.fun function
#
# Returns:
#  a list containing
#   - minimum coverage probability (aka. coverage coefficient)
#
######################################################################
coverage <- function(ci.fun, n, alpha=0.05, p.grid=NULL,interval=c(0,1),
                     pmfX=function(k,n,p) dbinom(k,size=n,prob=p), ...) {
  if (! (all(interval >= 0) & all(interval <= 1)) & (interval[2]<interval[1])) {
    stop("Error: Interval has to be a subset of [0,1] with borders upper>lower.")
  }

  #The minimum coverage probability is attained in the finite set of
  #probabilities that contains the left limit and right limit for
  #all \hat{p}_n = 0,1/n,\ldots, 1. 
  CIs <- ci.fun(x=0:n, n=n, conf.level=1-alpha,...)[,c("lower","upper")]

  #Indiciator function to investigate whether a specific p
  #is covered by the confidence interval or not
  Ikp <- function(k,p) {
    ifelse( p >= CIs[k+1,1] & p <= CIs[k+1,2], 1, 0)
  }
  #Compute coverage probability for a specific proportion p
  coverageprob.onep <- function(p) {
    k <- 0:n
    sum( Ikp(k,p) * pmfX(k,n,p) )  
  }

  #All values for p to evaluate
  p.grid <- c(p.grid,as.vector(as.matrix(CIs)))
  p.grid <- unique(sort(c(p.grid,seq(0,1,by=1/n))))
  p.grid <- subset(p.grid, p.grid >= interval[1] & p.grid <= interval[2])
  
  #Calculate coverage prob at each p
  probs <- sapply(p.grid, coverageprob.onep)

  #Coverage coefficient
  coeff <- min(probs)
  names(coeff) <- p.grid[which.min(probs)]

  #Done
  res <- list(coef=coeff,p=p.grid,probs=probs,alpha=alpha)
  class(res) <- "coverage"
  return(res)
}

plot.coverage <- function(x, y=NULL, ...) {
  plot(x$p, x$prob, xlab="p", ylab="coverage probability",...)
}

## ######################################################################
## # Find sample size, such that coverage coefficient in a certain
## # interval is just above a specific threshold
## #
## # This function only works at reasonable speed for small n.
## ######################################################################

## ciss.covcoeff <- function(ci.fun, alpha=0.05, cc.lowbound=(1-alpha)*0.95,limit=c(0,1),...) {

##   #Function to find root of. Cheat: n continous is rounded to give integer
##   f <- function(n,...) {
##     n <- round(n)
##     res <- coverage(ci.fun=binom.confint, n=n, alpha=alpha, interval=interval, ...)
##     cat("n= ",n,"\tCovCoef = ",res$coef,"\n")
##     diff <- res$coef - cc.lowbound
##     #Poor mans integer optimization
##     return( ifelse(diff>0, diff, 1e9*diff))
##   }

##   #Start at lowest value and advance
##   n <- 2
##   stop <- FALSE
##   while (!stop & (n < n0)) {
##     cat("n= ",n,"\tCovCoef = ",res$coef,"\n")
##     res <- coverage(ci.fun=binom.confint, n=n, alpha=alpha, method="logit",limit=c(1e-12,1-1e-12),interval=interval)
##     stop <- res$coef > cc.lowbound
##     n <- n+1
## #    plot(res$p,res$prob,type="s")
##   }
## }


#  optimres <- optimize(f, interval=c(100,1000),...)
#  n0 <- round(optimres$minimum)
#  all <- sapply(200:300, function(n) {
#    coverage(ci.fun=binom.confint, n=n, alpha=alpha, limit=limit,interval=interval,method="logit")$coef
#  })
#  f(100, method="logit")
#  plot(200:300,all - cc.lowbound)
#  
#  f(1000)
#  #Start value to begin optimization from
#  all <- sapply(2:n0, function(n) {
#    coverage(ci.fun=binom.confint, n=n, alpha=alpha, limit=limit,interval=interval,method="logit")$coef
#  })

#  lines(c(2,n0),rep(cc.lowbound,2))
#  plot(2:n0,all,type="l")
#  lines(c(2,n0),rep(cc.lowbound,2))
  
#cov <- coverage(ci.fun=binom.liubailey, n=256, alpha=0.1, lambda=5, d=0.05)
#cov <- coverage(ci.fun=binom.liubailey, n=256, alpha=0.1, p.grid=seq(0,1,length=1000),lambda=5, d=0.05)


#cov <- coverage(ci.fun=binom.confint, n=10, alpha=0.05, p.grid=seq(1e-12,1-1e-12,length=1000),method="wilson",interval=c(0.1,1))
#cov <- coverage(ci.fun=binom.confint, n=10, alpha=0.05, p.grid=seq(1e-12,1-1e-12,length=1000),method="logit")##

#cov <- coverage(ci.fun=binom.confint, n=10, alpha=0.05, p.grid=seq(1e-12,1-1e-12,length=1000),method="asymptotic")

#cov <- coverage(ci.fun=binom.liubailey, n=10, alpha=0.05, p.grid=seq(0,1,length=1000),lambda=1, d=0.1)
#cov$coef
#plot(cov$p, cov$probs,type="l")

#stop

## ##Function for constant width interval (strange interval!)
## Ikp.cw <- function(k, p, n, alpha, d) {
##   p.hat <- k/n
##   ci <- cbind(p.hat,p.hat) + outer(rep(d,length(k)),c(-1,1))
##   ifelse( p >= ci[,1] & p <= ci[,2], 1, 0)
## }

## ##Function for (3.1) intervals
## Ikp.liubailey <- function(k, p, n, alpha, lambda, d) {
##   ci <- binom.liubailey(k, n, lambda, d, conf.level = 1-alpha)[,c("lower","upper")]
##   ifelse( p >= ci[,1] & p <= ci[,2], 1, 0)
## }


## n <- 256
## cw <- cov.coeff(n=n, Ikp=Ikp.cw, d=0.05, alpha=0.1)
## cw$coef
## plot(cw$p.grid, cw$probs,type="l",xlab=expression(pi),ylab="coverage probability",ylim=c(0.7,1))
## lines(c(0,1),rep(1-cw$alpha,2),lty=2)


## liubailey <- cov.coeff(n=n, Ikp=Ikp.liubailey, d=0.05,lambda=5,alpha=0.1)
## liubailey$coef
## plot(liubailey$p.grid, liubailey$probs,type="l",xlab=expression(pi),ylab="coverage probability",ylim=c(0.7,1))
## lines(c(0,1),rep(1-liubailey$alpha,2),lty=2)


## d <- 0.05
## alpha <- 0.1
## (n0 <- ceiling( (0.5*qnorm(1-alpha/2)/d)^2))
## nstar <- n0
## stop <- FALSE
## while (!stop) {
##   nstar <- nstar + 1
##   print(nstar)
## #  liubailey <- cov.coeff(n=nstar, Ikp=Ikp.liubailey, d=d,lambda=5,alpha=alpha)
##   res <- cov.coeff(n=nstar, Ikp=Ikp.wilson, alpha=alpha)
##   print(res$coef)
##   stop <- res$coef >  1-alpha
## }
## #cov.prob(0.41,n=nstar,Ikp.cw,d=0.04)
## min(cov.liubailey)
## p.grid[which.min(cov.liubailey)]
