#******************************************************************************* 
#
# Estimation for Multivariate Normal Data with Monotone Missingness
# Copyright (C) 2007, University of Cambridge
# 
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
# 
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
# 
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
#
# Questions? Contact Robert B. Gramacy (bobby@statslab.cam.ac.uk)
#
#*******************************************************************************


## plot.blasso
##
## generic plotting method for blasso class objects.
## plots summaries of the samples from
## from the posterior distribution of the Bayesian lasso
## model

'plot.blasso' <-
  function(x, which=c("coef", "s2", "lambda2", "gamma",
                "tau2i", "omega2", "nu", "m", "pi"),
           subset=NULL, burnin=0, ...)
{
  ## check the burnin argument
  if(length(burnin) != 1 || burnin<0 || burnin >= x$T)
    stop("burnin must be a non-negative scalar < x$T")
  burnin <- burnin+1
  
  ## check the which argument
  which <- match.arg(which)

  ## make the appropriate kind of plot
  if(which == "coef") {
    if(is.null(subset) && x$M > 0) if(x$M > 0) subset <- 1:x$M
    if(!is.null(subset)) {
      beta <- x$beta[burnin:x$T,subset]
      df <- data.frame(mu=x$mu[burnin:x$T], beta)
    } else df <- data.frame(mu=x$mu[burnin:x$T])
    boxplot(df, ylab="coef", main="Boxplots of regression coefficients", ...)
    m <- which.max(x$lpost)
    points(c(x$mu[m], x$beta[m,]), col=2, pch=21)
    abline(h=0, lty=2, lwd=2)
  } else if(which == "s2") {
    par(mfrow=c(2,1))
    hist(x$s2[burnin:x$T], ...)
    plot(burnin:x$T, x$s2[burnin:x$T], type="l", main="s2 chain", ...)
  } else if(which == "lambda2") {
    if(is.null(x$lambda2)) stop("LASSO regression was not used")
    par(mfrow=c(2,1))
    hist(x$lambda2[burnin:x$T], ...)
    plot(burnin:x$T, x$lambda2[burnin:x$T], type="l", main="lambda2 chain", ...)
  } else if(which == "gamma") {
    if(is.null(x$gamma)) stop("NG regression was not used")
    par(mfrow=c(2,1))
    hist(x$gamma[burnin:x$T], ...)
    plot(burnin:x$T, x$gamma[burnin:x$T], type="l", main="gamma chain", ...)
  } else if(which == "m") {
    if(!x$RJ) stop("RJ was not used")
    par(mfrow=c(2,1))
    hist(x$m[burnin:x$T], ...)
    plot(burnin:x$T, x$m[burnin:x$T], type="l", main="m chain", ...)
  } else if(which == "tau2i"){
    if(is.null(subset)) subset <- 1:x$M
    if(is.null(x$tau2i)) stop("LASSO regression was not used")
    boxplot(data.frame(x$tau2i[burnin:x$T,subset]),
            main="Boxplot of tau2i", ylab="tau2i", ...)
  } else if(which == "omega2"){
    if(is.null(subset)) subset <- 1:x$n
    if(is.null(x$omega2)) stop("Student-t errors were not used")
    boxplot(data.frame(x$omega2[burnin:x$T,subset]),
            main="Boxplot of omega2", ylab="omega2", ...)
  } else if(which == "nu") {
    if(is.null(x$nu)) stop("Student-t errors were not used")
    par(mfrow=c(2,1))
    hist(x$nu[burnin:x$T], ...)
    plot(burnin:x$T, x$nu[burnin:x$T], type="l", main="nu chain", ...)
  } else if(which == "pi"){
    if(is.null(x$pi)) stop("the prior held pi fixed")
    par(mfrow=c(2,1))
    hist(x$pi[burnin:x$T], ...)
    plot(burnin:x$T, x$pi[burnin:x$T], type="l", main="pi chain", ...)
  }
}


## summary.blasso
##
## generic summary method for blasso class objects,
## basically calls summary on the matrices and vectos
## of samples from the posterior distribution of the
## parameters in the Bayesian lasso model

'summary.blasso' <-
function(object, burnin=0, ...)
  {

    ## check the burnin argument
    if(length(burnin) != 1 || burnin<0 || burnin >= object$T)
      stop("burnin must be a non-negative scalar < object$T")
    burnin <- burnin+1

    ## make the list
    rl <- list(call=object$call, B=burnin-1, T=object$T, thin=object$thin)
    class(rl) <- "summary.blasso"
    
    ## call summary on each object
    df <- data.frame(mu=object$mu[burnin:object$T],
                     object$beta[burnin:object$T,])
    rl$coef <- summary(df)
    rl$s2 <- summary(object$s2[burnin:object$T])

    ## only do if lasso, hs or ng
    if(!is.null(object$lambda2))
      rl$lambda2 <- summary(object$lambda2[burnin:object$T])
    if(!is.null(object$tau2i))
      rl$tau2i <- summary(data.frame(object$tau2i[burnin:object$T,]))

    ## only do if ng prior
    if(!is.null(object$gamma))
      rl$gamma <- summary(object$gamma[burnin:object$T])
    
    ## only do if Student-t errors
    if(!is.null(object$omega2))
      rl$omega2 <- summary(data.frame(object$omega2[burnin:object$T,]))
    
    ## only do if Student-t errors
    if(!is.null(object$nu))
      rl$nu <- summary(data.frame(object$nu[burnin:object$T]))
    
    ## only do if RJ
    if(object$RJ) {
      rl$bn0 <- apply(as.matrix(object$beta[burnin:object$T,]), 2,
                      function(x){ sum(x != 0) })/(object$T-burnin+1)
      rl$m <- summary(object$m[burnin:object$T])
    }

    ## only do if pi not fixed
    if(!is.null(object$pi))
      rl$pi = summary(object$pi[burnin:object$T])
    
    ## print it or return it
    rl
  }


## print.summary.blasso
##
## print the results of the summary method after first
## calling the print method on the blasso object.  

'print.summary.blasso' <-
  function(x, ...)
{
  ## print information about the call
  cat("\nCall:\n")
  print(x$call)
  
  ## print the monomvn object
  cat("\nsummary of MCMC samples with B=", x$B, " burnin rounds\n", sep="")
  cat("T=", x$T, " total rounds, with thin=", x$thin,
      " rounds between\n", sep="")
  cat("each sample\n\n")
  
  ## print coef
  cat("coefficients:\n")
  print(x$coef)
  cat("\n")

  ## print s2
  cat("s2:\n")
  print(x$s2)
  cat("\n")

  ## print lambda2
  if(!is.null(x$lambda2)) {
    cat("lambda2:\n")
    print(x$lambda2)
    cat("\n")
  }

  ## print gamma
  if(!is.null(x$gamma)) {
    cat("gamma:\n")
    print(x$gamma)
    cat("\n")
  }

  ## print tau2i
  if(!is.null(x$tau2i)) {
    cat("tau2i:\n")
    print(x$tau2i)
    cat("\n")
  }

  ## print bn0
  if(!is.null(x$bn0)) {
    cat("probability of beta != 0:\n")
    print(x$bn0)
    cat("\n")
  }

  ## print pi
  if(!is.null(x$pi)) {
    cat("pi:\n")
    print(x$bn0)
    cat("\n")
  }
}


## print.blasso
##
## generic print method for blasso class objects,
## summarizing the results of a blasso call

`print.blasso` <-
function(x, ...)
  {
    ## print information about the call
    cat("\nCall:\n")
    print(x$call)

    ## print the monomvn object
    cat("\nrun for T=", x$T, " MCMC samples, with thin=", x$thin,
        " rounds\nbetween each sample\n", sep="")

    ## say something about lasso
    if(!is.null(x$tau2i)) {
      if(x$hs)
        cat("\nHorseshoe was used to shrink regression coefficients\n")
      else if(!is.null(x$gamma))
        cat("\nA Normal-Gamma prior was used to shrink regression coefficients\n")
      else cat("\nLasso was used to shrink regression coefficients\n")
    } else if(!is.null(x$lambda2)) {
      cat("\nA ridge parameter was used to shrink regression coefficients\n")
    }

    ## say something about model selection
    if(x$RJ) {
      cat("\nReversible Jump (RJ) was used to average over\n")
      cat("subsets of columns in the design matrix, using a\n")

      ## add in info about the prior
      if(x$mprior[1] == 0)
        cat("uniform prior on m in {0,...,", x$M, "}\n", sep="")
      else if(length(x$mprior) == 1)
        cat("Bin(m|n=", x$M, ",p=", x$mprior, ") prior\n", sep="")
      else cat("Bin(m|n=", x$M, ",p) prior with p~Beta(",
                 x$mprior[1], ",", x$mprior[2], ")\n", sep="")
        
    }

    ## suggestion
    cat("\nTry summary.blasso and plot.blasso on this object\n\n")
  }
