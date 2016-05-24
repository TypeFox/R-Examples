# TBSSurvival package for R (http://www.R-project.org)
# Copyright (C) 2013 Adriano Polpo, Cassio de Campos, Debajyoti Sinha
#                    Jianchang Lin and Stuart Lipsitz.
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

######
############
##################
########################
###########################################################################
## This file contains auxiliar functions - NOT to be used directly by users
###########################################################################

## on attach, just print the version number of the package
.onAttach <- function(lib,pkg)
{
  packageStartupMessage("TBSSurvival 1.2 loaded\n")
}

##  Density for mixture of uniform-exponential. This density is use as prior for TBS model, has between (a,b)
##  an uniform distribution and after `b' the tail is exponential. the parameter `p' dictates the volume of the
##  uniform part and `1-p' the volume of the tail.
##  \item{x}{vector of quantiles.}
##  \item{p}{number `0 < p < 1'.}
##  \item{a}{parameter of uniform, default=0.}
##  \item{b}{parameter that difine the end of uniforme part.}
## \value{ `dunif.exp' gives the density. }
.dunif.exp <- function(x,a=0,b,p) {
  out <- rep(0,length(x))
  for (i in 1:length(x)) {
    if (x[i] <= b)
      out[i] <- p*dunif(x[i],a,b)
    else
      out[i] <- dexp(x[i],-log(1-p)/b)
  }
  return(out)
}

## computes the log posterior function. Some priors are fixed within the method, as described in the paper:
## Transform both sides: a parametric approach, Polpo et al.
.logpost <- function(par,time,delta,dist,x=NULL,mean=5,sd=5) {
  lambda <- par[1]
  xi <- par[2]
  if (lambda > 0 && xi > 0) {
    beta <- par[3:length(par)]
    ## check if the arguments are all ok
    aux <- .test.tbs(lambda,xi,beta,x,time=time,type="d")
    ## .test.tbs may eventually re-cast beta and x, so we update them here
#    beta <- aux$beta
#    x <- aux$x

    if (dist$name != "t") {
      ## for any dist not t-student, use a mixture of uniform-exponential with the flat prior
      ## between 0 and 3 for lambda with 0.8 of the mass, a mixture of uniform-exponential with
      ## flat prior for xi between 0 and 2 with 0.9 of the mass, and a normal distribution with
      ## mean and sd as given by the arguments for the beta
      out <- log(.dunif.exp(lambda,a=0.00001,b=3,p=0.8))+
             log(.dunif.exp(xi,a=0.00001,b=2,p=0.9))+
             sum(log(dnorm(beta,mean,sd)))+
             .lik.tbs(par=par,time=time,delta=delta,dist=dist,x=x)
    } else {
      ## in the case of the t-student, there is no prior for xi to be included, but the others
      ## are as before
      out <- log(.dunif.exp(lambda,a=0.00001,b=3,p=0.8))+
             sum(log(dnorm(beta,mean,sd)))+
             .lik.tbs(par=par,time=time,delta=delta,dist=dist,x=x)
    }
  } else {
    ## if lambda is not positive and xi is not positive, the result is -inf
    out <- log(0)
  }
  ## any numerical issue sets the output to -inf
  if (is.nan(out) || is.na(out))
    out <- log(0)
  return(out)
}

## this function returns the current amount of spent time by
## the current R process (user time + system time) in minutes
.gettime <- function() {
  t <- proc.time()
  ## note that t[1]+t[2] are the process time, not the real computer time,
  ## so we still take proper values even if the machine is running other stuff
  out <- (t[1]+t[2])/60
  names(out) <- NULL
  return(out)
}

## computes the likelihood for the TBS model. par is an array with the lambda in the 1st
## position, xi in the 2nd and beta over the other positions of the array (starting form beta0,
## which must always exist no matter the number of covariates). time and delta are the
## survival/reliability information (1-delta is the censoring indicator), dist is the error
## distribution of the TBS model (has to be one of the available ones), x gives the information
## about the covariates, and notinf flags whether to turn -inf into a very negative number, which
## in some cases is useful for numerical reasons.
.lik.tbs <- function(par,time,delta,dist=dist.error("norm"),x=NULL,notinf=FALSE)
{
  ## split the info from par
  lambda <- par[1]
  xi     <- par[2]
  ## at least one beta must exist, so length(par) >= 3
  beta   <- par[3:length(par)]
  out <- -Inf
  ## result is -inf unless all vars below are positive
  if ((xi > 0.0001) && (all(time > 0)) && (lambda > 0.0001))
  {
    ## check if the arguments are all ok
    aux <- .test.tbs(lambda,xi,beta,x,time=time,type="d")
    ## .test.tbs may eventually re-cast beta and x, so we update them here
    beta <- aux$beta
    x <- aux$x

    out <- 0
    ## if there are covariates, R requires us treat x as a matrix, but that is the only
    ## difference between the if and else statements here (which we write just to avoid copying x
    ## to a temporary matrix
    if (is.matrix(x)) {
      ## for the cases where delta is one (that is, no censoring), just use the density function
      if (any(delta == 1)) {
        d.aux <- .dtbs(time=time[delta==1],lambda=lambda,xi=xi,beta=beta,x=x[delta==1,],dist=dist)
        out <- out + sum(log(d.aux))
      }
      ## for those where delta is zero (censored), use the survival function
      if (any(delta == 0)) {
        s.aux <- 1-.ptbs(time=time[delta==0],lambda=lambda,xi=xi,beta=beta,x=x[delta==0,],dist=dist)
        out <- out + sum(log(s.aux))
      }
    } else {
      ## as said before, this is a repetition of the "if" just because here x is not a matrix
      if (any(delta == 1)) {
        d.aux <- .dtbs(time=time[delta==1],lambda=lambda,xi=xi,beta=beta,x=x[delta==1],dist=dist)
        out <- out + sum(log(d.aux))
      }        
      if (any(delta == 0)) {
        s.aux <- 1-.ptbs(time=time[delta==0],lambda=lambda,xi=xi,beta=beta,x=x[delta==0],dist=dist)
        out <- out + sum(log(s.aux))
      }
    }
  }
  ## if it is not to return -inf, then return a very negative value... as we are dealing with logs, -1e10 suffices
  if(is.na(out) || length(out)==0) return(NA)
  if(notinf && out < -1e10) out = -1e10
  if(notinf && out > 1e10) out = 1e10
  return(out)
}

## nstart is the number of (feasible!) initial points to use. The method will try many guesses to find feasible points
## max.time (in minutes) to run the optimization
## method has to be one of the available in the function optim or "Rsolnp"
## dist has to be one of those available in the dist.error function (see file tbs.r)
## NOTICE: this function uses withTimeout from the R.utils package. We have experienced some versions of R.utils
##         which do not have this function (e.g. some versions installed with apt-get in ubuntu). In this case,
##         one has to install the CRAN version of R.utils
.tbs.survreg <- function(formula,dist=dist.error("norm"),method="BFGS",guess=NULL,nstart=10,verbose=FALSE,max.time=-1,gradient=TRUE) {
  initial.time <- .gettime()
  if(max.time <= 0) {
    ## user didn't define a timeout, so we set to a large number
    max.time <- 1e10
  }
  if (attributes(formula)$class != "formula")
    stop("A formula argument is required")
  grad=NULL
  if(gradient && (method != "SANN")) grad=.grad.tbs

  ## read the information from within the formula to populate the required variables
  mf <- model.frame(formula=formula)
  x <- model.matrix(attr(mf, "terms"), data=mf)
  y <- model.response(mf)
  time <- y[,1]
  delta <- y[,2]
  x.k   <- dim(x)[2]
  n     <- dim(x)[1]
  ## checks if delta is an indicator function
  if (any((delta != 0) & (delta != 1)))  {
    stop("Only uncesored or right censored data are allowed")
  }
  out <- NULL

  ## check if starting point was given or not, and build one in case. But note that
  ## we can't be sure this point will evaluate to something > -inf, so there might be
  ## the need later to search for more starting points...
  nparam <- 2
  if (!is.null(x)) {
    if (is.matrix(x))
      nparam <- nparam+length(x[1,])
    else
      nparam <- nparam+1
  }
  if(is.null(guess)) {
    ## betas can be anything, we sample uniformly from -10 to 10 (note this is arbitrary)
    guess <- 20*runif(nparam)-10
    ## lambda and xi have to be positive, and "desirable" values are not very high...
    guess[1] <- 5*runif(1)+0.0001 ## lambda (note 5 here is completely arbitrary)
    guess[2] <- 10*runif(1)+0.0001 ## xi (note 10 here is completely arbitrary)
  }
  if(nparam != length(guess))
    stop("Number of parameters in the formula and length of the initial guess do not match")

  if(method=="Rsolnp") {
      out$method <- method
      ## Rsolnp needs some bounds for the unknowns. We arbitrarily define them to be between -100 and 100.
      LB = rep(-100,nparam)
      UB = rep(100,nparam)
      ## unless for lambda and xi, which must always be positive
      LB[1] = 0.0001
      LB[2] = 0.0001
      ## upper bound for xi is not very clear, but 1000 should be enough. Lambda is left with 100 as upper.
      UB[2] = 1000
      if(verbose) cat('RSOLNP: ')
      for(itk in 1:3) {
        ans = try(withTimeout(gosolnp(pars = NULL, fixed = NULL,
          fun = function(pars, n) { -.lik.tbs(pars,time=time,delta=delta,x=x,dist=dist,notinf=TRUE) },
        LB = LB, UB = UB, control = list(outer.iter = 200, trace = 0, tol=1e-4, delta=1e-6),
        distr = rep(1, length(LB)), distr.opt = list(), n.restarts = nstart, n.sim=3000, rseed = runif(n=1,min=1,max=10000000), n = nparam),
        timeout=max.time*60,onTimeout="error"))
        if (class(ans) != "try-error" && ans$convergence == 0 && length(ans$values)>0 && ans$values[length(ans$values)] < 1e10) {
          break
        }
######        print(ans)
      }
      if (class(ans) != "try-error" &&  ans$convergence == 0 && length(ans$values)>0 && ans$values[length(ans$values)] < 1e10) {
        ## process the solution in case one was found
        ## get parameters
#        out$par <- ans$pars
        out$lambda <- ans$par[1]
        out$xi <- ans$par[2]
        out$beta <- ans$par[3:length(ans$par)]
        options(warn = -1)
        ## compute the std.error
        aux <- try(sqrt(diag(solve(-(ans$hessian)))),silent=TRUE)
        options("warn" = 0)
        std.error <- rep(NA,nparam)
        if (class(aux) != "try-error")
          std.error <- aux
        out$lambda.se <- std.error[1]
        out$xi.se <- std.error[2]
        out$beta.se <- std.error[3:length(std.error)]
        ## get the log-lik value
        out$log.lik <- -ans$values[length(ans$values)]
        if(verbose) cat(out$log.lik,'PARS:',ans$pars,'TIME:',ans$elapsed,'\n')
        out$error.dist <- dist
        ## compute error distances
        out$AIC  <- 2*nparam-2*out$log.lik
        out$AICc <- 2*nparam-2*out$log.lik + 2*nparam*(nparam+1)/(length(time)-nparam-1)
        out$BIC  <- -2*out$log.lik+nparam*log(length(time))
        out$convergence <- TRUE
        ## evaluate the "error"
        aux <- .test.tbs(out$lambda,out$xi,out$beta,x,time,type="d")
#        out$time  <- time[delta == 1]
#        out$error <- c(.g.lambda(log(out$time),out$lambda)-.g.lambda(c(aux$x%*%aux$beta)[delta == 1],out$lambda))
        out$time  <- time
        out$delta <- delta[order(time)]
        out$error <- c(.g.lambda(log(out$time),out$lambda)-.g.lambda(c(aux$x%*%aux$beta),out$lambda))[order(time)]
        out$time  <- time[order(time)]

        names(out$time) <- NULL
        names(out$error) <- NULL
        ## for the plot
        if (length(out$beta) == 1) {
          if (unique(aux$x[,1]) == 1) {
            out$x <- 1
            attr(out$x,"plot") <- 1
          } else if (length(unique(aux$x[,1]) <= 4)) {
            out$x <- unique(aux$x[,1])
            attr(out$x,"plot") <- 2
          }
        } else if ((length(out$beta) == 2) && (unique(aux$x[,1]) == 1) &&
                   (length(unique(aux$x[,2])) <= 4)) {
          out$x <- unique(aux$x[,2])
          attr(out$x,"plot") <- 3
        }
        ## set time run time
        out$run.time <- .gettime() - initial.time
      } else {
        if(verbose) cat(' failed\n')
        out$run.time <- .gettime() - initial.time
        out$convergence <- FALSE
        cat(paste(method,": It was not possible to find a feasible solution\n"))
      }
      return(out)
  }

  ## i will count the number of feasible guesses already used
  i <- 1
  ## est will keep the best solution found so far
  est=NA
  ## ii will count the number of unfeasible initial guesses
  ii=1
  if(verbose) cat(method,': ',sep='')
  inimethod=method
  wasnan=TRUE
  ## we also control the maximum amount of time this method can keep looking for better solutions
  inilooptime=.gettime()
  while(.gettime() < inilooptime + max.time) {
      valik=.lik.tbs(guess,time=time,delta=delta,x=x,dist=dist)
    ## check if the guess evaluates to -inf, in this case it is not worth to spend time in the optim, unless
    ## we have not found any feasible point yet. In this case, better try it...
    if(!is.na(valik) && (valik>-Inf || is.na(est))) {
      aux <- try(withTimeout(optim(guess, fn=.lik.tbs, gr=grad, time=time, delta=delta, dist=dist, x=x, notinf=FALSE,
                                   method=inimethod, control=list(fnscale=-1), hessian=TRUE),timeout=max.time*60,onTimeout="error"),silent=TRUE)
      if (class(aux) != "try-error") {
##              print(paste('aux$value',aux$value))
        if((inimethod=="SANN") || (aux$convergence != 0)) {
          for(itx in 1:10) {
            aux1 <- try(withTimeout(optim(aux$par, fn=.lik.tbs, gr=grad, time=time, delta=delta, dist=dist, x=x, notinf=FALSE,
                                          method=method, control=list(fnscale=-1), hessian=TRUE),timeout=max.time*60,onTimeout="error"),silent=TRUE)
##                print(paste('aux1$value',aux1$value))
            if (class(aux1) != "try-error") {
              if (aux1$value < aux$value + 0.0001) {
                ## 0.0001 is only for numerical reasons. Note that 0.0001 in the log value is anyway very very small...
                break
              }
              aux=aux1
            } else {
              break
            }
          }
        }
        ## if a new best solution was found, update the current one
        if(is.na(est) || aux$value > est$value || wasnan) {
          options("warn" = -1)
          aux1 <- try(sqrt(diag(solve(-(aux$hessian)))),silent=TRUE)
          options("warn" = 0)
          std.error <- rep(NaN,nparam)
          if (class(aux1) != "try-error")
            std.error <- aux1
          willnan <- is.nan(sum(std.error))
          if(is.na(est) || (wasnan && !willnan) || (aux$value > est$value && wasnan==willnan)) {
            wasnan = is.nan(sum(std.error))
            est = aux
            inimethod=method
          }
        }
        ## one more feasible point found
        i = i + 1
        if(verbose) cat('@')
      } else {
##        print(aux)
        if(verbose) cat('*')
        ## one more unfeasible point tried
        ii = ii + 1
      }
    } else {
      ## one more unfeasible point tried
      ii = ii + 1
    }
    ## betas can be anything, we sample uniformly from -10 to 10
    guess <- 20*runif(nparam)-10
    ## lambda and xi have to be positive, and "desirable" values are not very high...
    guess[1] <- 5*runif(1)+0.0001 ## lambda
    guess[2] <- 10*runif(1)+0.0001 ## xi

    ## if enough starts have been tried, stop. Also stop if too many unsuccessfull tries have been made :(
    if(ii>100 && is.na(est)) {
      inimethod="SANN"
      grad=NULL
      if(verbose) cat('$')
    }
    ## ii counts the number of initial guess that have been tried but were unfeasible ones. We try at least
    ## one thousand times even if the users selected fewer nstart, because nstart is usually meant to be the
    ## number of feasible initial guess...
    if(i>nstart || ii > max(nstart,1000)) {
      break
    }
  }

  out$method <- method
  if(!is.na(est) && est$value > -Inf) {
    ## if there is at least one feasible solution that has been found, compute some quantities for it,
    ## such as AIC, BIC, std.err, and return them
#    out$par <- est$par
    out$lambda <- est$par[1]
    out$xi <- est$par[2]
    out$beta <- est$par[3:length(est$par)]
    options(warn = -1)
    ## compute the std.error
    aux <- try(sqrt(diag(solve(-(est$hessian)))),silent=TRUE)
    options("warn" = 0)
    std.error <- rep(NA,nparam)
    if (class(aux) != "try-error")
      std.error <- aux
    out$lambda.se <- std.error[1]
    out$xi.se <- std.error[2]
    out$beta.se <- std.error[3:length(std.error)]
    if(is.nan(sum(std.error))) {
      warning(paste('tbs.survreg.mle: optimization method',method,'failed to compute standard errors -- possibly it is not in a local optimum'))
    }
    ## get the log.lik value
    out$log.lik <- est$value
    out$error.dist <- dist
    ## evaluate some scores
    out$AIC  <- 2*nparam-2*est$value
    out$AICc <- 2*nparam-2*est$value + 2*nparam*(nparam+1)/(length(time)-nparam-1)
    out$BIC  <- -2*est$value+nparam*log(length(time))
    out$convergence <- TRUE
    ## evaluate the "error"
    aux <- .test.tbs(out$lambda,out$xi,out$beta,x,time,type="d")
##    out$time  <- time[delta == 1]
##    out$error <- c(.g.lambda(log(out$time),out$lambda)-.g.lambda(c(aux$x%*%aux$beta)[delta == 1],out$lambda))
    out$time  <- time
    out$delta <- delta[order(time)]
    out$error <- c(.g.lambda(log(out$time),out$lambda)-.g.lambda(c(aux$x%*%aux$beta),out$lambda))[order(time)]
    out$time  <- time[order(time)]
    names(out$time) <- NULL
    names(out$error) <- NULL
    ## for the plot
    if (length(out$beta) == 1) {
      if (unique(aux$x[,1]) == 1) {
        out$x <- 1
        attr(out$x,"plot") <- 1
      } else if (length(unique(aux$x[,1]) <= 4)) {
        out$x <- unique(aux$x[,1])
        attr(out$x,"plot") <- 2
      }
    } else if ((length(out$beta) == 2) && (unique(aux$x[,1]) == 1) &&
               (length(unique(aux$x[,2])) <= 4)) {
      out$x <- unique(aux$x[,2])
      attr(out$x,"plot") <- 3
    }
    ## record run time
    out$run.time <- .gettime() - initial.time
    if(verbose) cat(' ',out$log.lik,'TIME:',out$run.time,'\n')
  } else {
    if(verbose) cat(' failed\n')
    out$convergence <- FALSE
    out$run.time <- .gettime() - initial.time
    cat(paste(method,": It was not possible to find a feasible solution\n"))
  }
  return(out)
}

.mylog <- function(t) log(t+1e-8)
.invn <- function(t) (1.0/(t+1e-8))
.pown <- function(a,b) exp(b*log(abs(a)+1e-8))

## \code{.g.lambda} gives the generalized power transformation function
.g.lambda <- function(x,lambda) {
  return(sign(x)*exp(lambda*log(abs(x)))/lambda)
#  return(sign(x)*(.pown(x,lambda)* .invn(lambda)))
}

## \code{.g.lambda.inv} is the inverse of generalized power transformation function.
.g.lambda.inv <- function(x,lambda) {
  return(sign(x)*(exp(log(abs(x*lambda))/lambda)))
#  return(sign(x)*(.pown(x*lambda,.invn(lambda))))
}

.deriv.tbs <- function(x, xi, dist, type, var) {
  switch(dist$name,
    norm = switch(var,
      xi = switch(type,
        dens = ((x^2)/(xi^(5/2))-1/(xi^(3/2)))*exp(-(x^2)/(2*xi))/sqrt((2^3)*pi),
        surv = -x*exp(-(x^2)/(2*xi))/sqrt(((2*xi)^3)*pi)),
      x  = switch(type,
        dens = -x*exp(-(x^2)/(2*xi))/sqrt(2*pi*xi^3),
        surv = -dist$d(x,xi))),
   t = switch(var,
      xi = switch(type,
        dens = ((((x^2)/xi+1)^((-xi-1)/2)*
                gamma((xi+1)/2)/(sqrt(pi*xi)*gamma(xi/2)))*
                (digamma((xi+1)/2)/2-digamma((xi)/2)/2
                 -(x^2)*(-xi-1)/(2*((x^2)/xi+1)*xi^2)-log((x^2)/xi+1)/2-1/(2*xi))),
        surv = ((-digamma((xi+1)/2)+digamma(xi/2)+(1/xi))*
                f21hyper(1/2,(xi+1)/2,3/2,-(x^2)/xi)*x*gamma((xi+1)/2)/(2*sqrt(pi*xi)*gamma(xi/2))
                -f21hyper(3/2,(xi+1)/2+1,5/2,-(x^2)/xi)*(x^3)*(xi+1)*gamma((xi+1)/2)/(6*sqrt(pi*xi^5)*gamma(xi/2)))),
      x  = switch(type,
        dens = x*(((x^2)/xi+1)^((-xi-1)/2-1))*(-xi-1)*gamma((xi+1)/2)/(sqrt(pi*xi^3)*gamma(xi/2)),
        surv = -dist$d(x,xi))),
    cauchy = switch(var,
      xi = switch(type,
        dens = ((2*x^2)/(((x^2)/(xi^2)+1)*(xi^2))-1)/(pi*((x^2)/(xi^2)+1)*(xi^2)),
        surv = x/(pi*((x^2)/(xi^2)+1)*(xi^2))),
      x  = switch(type,
        dens = -2*x/(pi*(((x^2)/(xi^2)+1)^2)*xi^3),
        surv = -dist$d(x,xi))),
    doubexp = switch(var,
      xi = switch(type,
        dens = (abs(x)/xi-1)*exp(-abs(x)/xi)/(2*xi^2),
        surv = (x/2)*exp(-abs(x)/xi)/(xi^2)),
      x  = switch(type,
        dens = -x*exp(-abs(x)/xi)/(2*abs(x)*xi^2),
        surv = -dist$d(x,xi))),
    logistic = switch(var,
      xi = switch(type,
        dens = (2*x*exp(2*x/xi)/((xi^3)*(exp(x/xi)+1)^3)-exp(x/xi)/((xi^2)*(exp(x/xi)+1)^2)
               -x*exp(x/xi)/((xi^3)*(exp(x/xi)+1)^2)),
        surv = -x*exp(-x/xi)/((xi^2)*(exp(-x/xi)+1)^2)),
      x  = switch(type,
        dens = exp(x/xi)/((xi^2)*(exp(x/xi)+1)^2)-2*exp(2*x/xi)/((xi^2)*(exp(x/xi)+1)^3),
        surv = -dist$d(x,xi))))
}

.deriv.logL.tbs <- function(t, eta, lambda, xi, dist, type, var) {
  switch(type,
    Se = switch(var,
      eta    = (- .pown(eta,lambda-1))*
                dist$d(.g.lambda(.mylog(t),lambda)-.g.lambda(eta,lambda),xi)*
                (.invn(1-dist$p(.g.lambda(.mylog(t),lambda)-.g.lambda(eta,lambda),xi))),
      lambda = (((sign(.mylog(t))* .pown(.mylog(t),lambda)/lambda)*(.mylog(abs(.mylog(t)))-1/lambda)
                -(sign(eta)* .pown(eta,lambda)/lambda)*(.mylog(abs(eta))-1/lambda))*
                (-dist$d(.g.lambda(.mylog(t),lambda)-.g.lambda(eta,lambda),xi))*
               (.invn(1-dist$p(.g.lambda(.mylog(t),lambda)-.g.lambda(eta,lambda),xi)))),
      xi     =  .deriv.tbs(.g.lambda(.mylog(t),lambda)-.g.lambda(eta,lambda),xi,dist,type="surv",var="xi")* 
                (.invn(1-dist$p(.g.lambda(.mylog(t),lambda)-.g.lambda(eta,lambda),xi)))),
    fe = switch(var,
      eta    = .pown(eta,lambda-1)*
                .deriv.tbs(.g.lambda(.mylog(t),lambda)-.g.lambda(eta,lambda),xi,dist,type="dens",var="x")*
                  (.invn(dist$d(.g.lambda(.mylog(t),lambda)-.g.lambda(eta,lambda),xi))),
      lambda = (((sign(.mylog(t))* .pown(.mylog(t),lambda)/lambda)*(.mylog(abs(.mylog(t)))-1/lambda)
                -(sign(eta)* .pown(eta,lambda)/lambda)*(.mylog(abs(eta))-1/lambda))*
                (.deriv.tbs(.g.lambda(.mylog(t),lambda)-.g.lambda(eta,lambda),xi,dist,type="dens",var="x"))*
                   (.invn(dist$d(.g.lambda(.mylog(t),lambda)-.g.lambda(eta,lambda),xi)))),
      xi     =  .deriv.tbs(.g.lambda(.mylog(t),lambda)-.g.lambda(eta,lambda),xi,dist,type="dens",var="xi")*
                  (.invn(dist$d(.g.lambda(.mylog(t),lambda)-.g.lambda(eta,lambda),xi)))))
}

.grad.tbs <- function(par,time,delta,dist,x=NULL,notinf=TRUE) {
  lambda <- par[1]
  xi     <- par[2]
  if(lambda < 0 || xi < 0) {
    if(is.matrix(x)) out <- rep(-Inf,2+length(x[1,]))
    else out <- c(-Inf,-Inf,-Inf)
    return(out)
  }
  ## at least one beta must exist, so length(par) >= 3
  beta   <- par[3:length(par)]

  if (is.null(x)) {
    x <- rep(1,length(time))
    eta <- beta*x
  } else {
    if (is.matrix(x)) {
      eta <- c(beta%*%t(x))
    } else {
      eta <- beta*x
    }
  }
  ## print(time)
  ## print(delta)
  ## print(par)
  ## print(x)
  
  aux.eta    <- rep(NA,length(time))
  aux.lambda <- rep(NA,length(time))
  aux.xi     <- rep(NA,length(time))
  for (i in 1:length(time)) {
    aux.eta[i] <-    (delta[i]*.deriv.logL.tbs(time[i],eta[i],lambda,xi,dist,type="fe",var="eta")+
                  (1-delta[i])*.deriv.logL.tbs(time[i],eta[i],lambda,xi,dist,type="Se",var="eta"))
    aux.lambda[i] <- (delta[i]*.deriv.logL.tbs(time[i],eta[i],lambda,xi,dist,type="fe",var="lambda")+
                  (1-delta[i])*.deriv.logL.tbs(time[i],eta[i],lambda,xi,dist,type="Se",var="lambda"))
    aux.xi[i] <-     (delta[i]*.deriv.logL.tbs(time[i],eta[i],lambda,xi,dist,type="fe",var="xi")+
                  (1-delta[i])*.deriv.logL.tbs(time[i],eta[i],lambda,xi,dist,type="Se",var="xi"))
  }
  ## print(aux.eta)
  ## print(aux.lambda)
  ## print(aux.xi)
  
  if (is.matrix(x)) {
    out <- rep(NA,2+length(x[1,]))
    out[1] <- sum(aux.lambda)
    out[2] <- sum(aux.xi)
    for (i in 1:length(x[1,]))
      out[i+2] <- sum(x[,i]*aux.eta)
  } else {
    out <- c(sum(aux.lambda),sum(aux.xi),sum(x*aux.eta))
  }
  if(notinf) {
    for(i in 1:length(out)) {
      if(out[i] < -1e10) out[i] = -1e10
      if(out[i] > 1e10) out[i] = 1e10
    }
  }
  return(out)
}

## this function has the sole purpose of checking whether the arguments respect the
## needs of the other TBS functions' implementation. It also re-cast the arguments in
## case it is needed, but does not really perform calculations.
.test.tbs <- function(lambda, xi, beta, x=NULL, time=NULL, type=NULL, p=NULL, n=NULL) {
  if (!is.numeric(xi))
    stop("xi is not a number")
  if (is.matrix(xi))
    stop("xi is matrix")
  if (xi <= 0)
    stop("xi <= 0")
  
  out   <- NULL
  out$x <- x
  out$beta <- beta

  if ((!is.numeric(lambda)) || (length(lambda) != 1))
    stop("lambda is not a number or length != 1")
  if (!is.numeric(beta))
    stop("beta is not a (vector) number")
  if (is.matrix(beta))
    stop("beta is matrix")
  if (!is.null(x)) {
    if (is.matrix(x)) {
      if (length(beta) != length(x[1,]))
        stop(paste("size of beta != ",length(x[1,]),sep=""))
    }
    else {
      if ((length(beta) != 1) && (length(beta) != length(x)))
        stop("size of beta is not conform")
    }
  }
  else {
    if (length(beta) > 1)
      stop("x is wrong or length(beta) > 1")  
  }
  if (lambda <= 0)
    stop("lambda <= 0")

  if (!is.null(type)) {
    if ((type == "d") || (type == "p")) {
      if (!is.numeric(time))
        stop("time is not a (vector) number")
      if (is.matrix(time))
        stop("time is matrix")
      if (any(time <= 0))
        stop("time <= 0")
      if (!is.null(x)) {
        if (is.matrix(x)) {
          if (length(time) != length(x[,1]))
            stop("length of time is different of length of x")
        }
        else {
          if (length(beta) == length(x)) {
            out$x <- matrix(x,1,length(x))
          } else {
            if (length(time) != length(x))
              stop("length of time is different of length of x")
            out$x <- matrix(x,length(x),1)
          }
        }
        out$beta <- matrix(beta,length(beta),1)
      }
      else {
        out$x <- matrix(1,length(time),1)
      }
    } else {
      if (type == "q") {
        if (!is.numeric(p))
          stop("p is not a (vector) number")
        if (is.matrix(p))
          stop("p is matrix")
        if (min(p) < 0)
          stop("p < 0")
        if (max(p) > 1)
          stop("p > 1")
      } else if (type == "r") {
          if (!is.numeric(n))
            stop("n is not a number")
          if (n %% 1 != 0)
            stop("n is not a integer number")
        }
      if (is.null(x)) {
        if (length(beta) > 1)
          stop("If x is omitted then beta must have length 1")
        out$x <- 1
      }
    }
  }

  return(out)
}
