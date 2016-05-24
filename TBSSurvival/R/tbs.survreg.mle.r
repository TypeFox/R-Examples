## TBSSurvival package for R (http://www.R-project.org)
## Copyright (C) 2013 Adriano Polpo, Cassio de Campos, Debajyoti Sinha
##                    Jianchang Lin and Stuart Lipsitz.
##
##    This program is free software: you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation, either version 3 of the License, or
##    (at your option) any later version.
##
##    This program is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.
##
##    You should have received a copy of the GNU General Public License
##    along with this program.  If not, see <http://www.gnu.org/licenses/>.

## Maximum likelihood estimation for TBS
## By default, try all optimization methods listed below and keep the best solution
## max.time is the time limit in minutes for each method (<= 0 means no limit), nstart the number of feasible starting points to use,
## dist defines the error distribution by name ("norm", "doubexp", "t", "cauchy", "logistic") or using
## a call from dist.error, or even created directly by the user (see dist.error.r for details)
## formula is a R formula with a Surv object on the left side
tbs.survreg.mle <- function(formula,dist=dist.error("all"),method=c("Nelder-Mead","BFGS","Rsolnp","SANN","CG"),
                            verbose=FALSE,nstart=10,max.time=-1,seed=1234,gradient=FALSE) {
  if(is.character(dist)) dist <- dist.error(dist)
  if(length(method)==1 && method=='all') method=c("Nelder-Mead","BFGS","Rsolnp","SANN","CG")
  ## this meta-method only records the elapsed time and call the max likelihood estimation function
  ## for each of the methods given until one of them converges. It is supposed that at least one method
  ## is given, and that dist is one of those implemented by tbs.survreg.
  initial.time <- .gettime()
  set.seed(seed)

  fn.aux <- function(formula,dist,method,verbose,nstart,max.time) {
    initial.time <- .gettime()
    bestout <- NULL
    for(i in 1:length(method)) {
      ## call the estimation function
      out <- .tbs.survreg(formula,dist=dist,method=method[i],verbose=verbose,max.time=max.time,nstart=nstart,gradient=gradient)
      ## if converged, we are happy
      if(out$convergence) {
        if(is.null(bestout)) {
          bestout <- out
        } else {
          wasnan <- is.nan(bestout$lambda.se) || is.nan(bestout$xi.se) || is.nan(sum(bestout$beta.se))
          isnan <- is.nan(out$lambda.se) || is.nan(out$xi.se) || is.nan(sum(out$beta.se))
          if (out$log.lik > bestout$log.lik && (wasnan || !isnan)) {
            bestout <- out
          }
        }
      }
    }
    ## check if at least one method found a solution...
    if (is.null(bestout)) {
      if(verbose) cat('\nTBS: No method has found a solution\n')
      bestout$convergence <- FALSE
      bestout$method <- NULL
    }
    ## record the call arguments and the formula
    bestout$call <- match.call()
    bestout$formula <- formula
    ## run.time is returned in minutes
    bestout$run.time <- .gettime() - initial.time
    
    class(bestout) <- "tbs.survreg.mle"
    return(bestout)
  }
  if (is.list(dist) && typeof(dist[[1]])=="list") {
    out <- NULL
    aux <- Inf
    for (i in 1:length(dist)) {
      if (verbose) cat(dist[[i]]$name,":\n",sep="")
      eval(parse(text=paste("out$",dist[[i]]$name," <- 1",sep=""))) 
      out[[i]] <- fn.aux(formula=formula,dist[[i]],method,verbose,nstart,max.time)
      if ((aux > out[[i]]$AIC) && (out[[i]]$convergence)) {
        aux <- out[[i]]$AIC
        best <- dist[[i]]$name
        bestn <- i
      }
    }
    if (exists("best")) {
      out$best <- best
      out$best.n <- bestn
      class(out) <- "tbs.survreg.mle.best"
    } else {
      out$best <- "none"
      out$best.n <- 0
    }
  } else {
    out <- fn.aux(formula,dist,method,verbose,nstart,max.time)
  }
  
  out$run.time <- .gettime() - initial.time
  return(out)
}

print.tbs.survreg.mle.best <- function(x, ...) {
  print(unclass(x), ...)
}

print.tbs.survreg.mle <- function(x, ...) {
#  print(x$call)
  if (x$convergence) {
    if (x$error.dist$name == "norm")
      text.dist <- "normal"
    if (x$error.dist$name == "t")
      text.dist <- "t-student"
    if (x$error.dist$name == "cauchy")
      text.dist <- "Cauchy"
    if (x$error.dist$name == "doubexp")
      text.dist <- "Double exponential"
    if (x$error.dist$name == "logistic")
      text.dist <- "logistic"
    else text.dist <- x$error.dist$name 
    cat("\n",sep="")
    cat("--------------------------------------------------------\n",sep="")
    cat(" TBS model with ",text.dist," error distribution (MLE).\n",sep="")
    cat("\n",sep="")
    
    aux1 <- nchar(sprintf("%.4f",c(x$lambda,x$xi,x$beta)))
    aux2 <- nchar(sprintf("%.4f",c(x$lambda.se,x$xi.se,x$beta.se)))
    text <- c("Estimates","Std. Error")
    auxc <- nchar(text)
    auxc1 <- c(ifelse(max(aux1) > auxc[1],abs(max(aux1)-auxc[1]),0),
               ifelse(max(aux2) > auxc[2],abs(max(aux2)-auxc[2]),0))
    auxc <- auxc+auxc1
    
    cat("         ",rep(" ",auxc1[1]+1),text[1],
        rep(" ",auxc1[2]+1),text[2],"\n",sep="")
    cat("  lambda: ",
        rep(ifelse(auxc[1] > aux1[1]," ",""),abs(auxc[1]-aux1[1])),sprintf("%.4f",x$lambda)," ",
        rep(ifelse(auxc[2] > aux2[1]," ",""),abs(auxc[2]-aux2[1])),sprintf("%.4f",x$lambda.se),
        "\n",sep="")
    cat("      xi: ",
        rep(ifelse(auxc[1] > aux1[2]," ",""),abs(auxc[1]-aux1[2])),sprintf("%.4f",x$xi)," ",
        rep(ifelse(auxc[2] > aux2[2]," ",""),abs(auxc[2]-aux2[2])),sprintf("%.4f",x$xi.se),
        "\n",sep="")
    if (length(x$beta) == 1) {
      cat("    beta: ",
          rep(ifelse(auxc[1] > aux1[3]," ",""),abs(auxc[1]-aux1[3])),sprintf("%.4f",x$beta)," ",
          rep(ifelse(auxc[2] > aux2[3]," ",""),abs(auxc[2]-aux2[3])),sprintf("%.4f",x$beta.se)," ",
          "\n",sep="")
    } else {
      for (i in 1:length(x$beta)) {
        if (i-1 < 10) {
          cat("   beta",i-1,": ",
              rep(ifelse(auxc[1] > aux1[(i+2)]," ",""),abs(auxc[1]-aux1[(i+2)])),sprintf("%.4f",x$beta[i])," ",
              rep(ifelse(auxc[2] > aux2[(i+2)]," ",""),abs(auxc[2]-aux2[(i+2)])),sprintf("%.4f",x$beta.se[i])," ",
              "\n",sep="")
        } else {
          cat("  beta",i-1,": ",
              rep(ifelse(auxc[1] > aux1[(i+2)]," ",""),abs(auxc[1]-aux1[(i+2)])),sprintf("%.4f",x$beta[i])," ",
              rep(ifelse(auxc[2] > aux2[(i+2)]," ",""),abs(auxc[2]-aux2[(i+2)])),sprintf("%.4f",x$beta.se[i])," ",
              "\n",sep="")
        }
      }
    }
    cat("\n",sep="")
    cat("     AIC: ",sprintf("%.4f",x$AIC),"\n",sep="")
    cat("    AICc: ",sprintf("%.4f",x$AICc),"\n",sep="")
    cat("     BIC: ",sprintf("%.4f",x$BIC),"\n",sep="")
    cat("\n",sep="")
    cat("Run time: ", sprintf("%.2f",x$run.time)," min \n",sep="")
    cat("--------------------------------------------------------\n",sep="")
    cat("\n",sep="")
  } else {
    cat("There was no convergence.\n",sep="")
  }
}

summary.tbs.survreg.mle.best <- function(object, ...) {
  summary.tbs.survreg.mle(object[[object$best.n]], ...)
}

summary.tbs.survreg.mle <- function(object, ...) {
  x=object
  if (x$convergence) {
    if (x$error.dist$name == "norm")
      text.dist <- "normal"
    if (x$error.dist$name == "t")
      text.dist <- "t-student"
    if (x$error.dist$name == "cauchy")
      text.dist <- "Cauchy"
    if (x$error.dist$name == "doubexp")
      text.dist <- "Double exponential"
    if (x$error.dist$name == "logistic")
      text.dist <- "logistic"
    else text.dist <- x$error.dist$name 
    cat("--------------------------------------------------------\n",sep="")
    cat(" TBS model with ",text.dist," error distribution (MLE).\n",sep="")
    cat("\n",sep="")
    
    z.value <- rep(NA,length(x$beta))
    p.value <- rep(NA,length(x$beta))
    for (i in 1:length(x$beta)) {
      if ((!is.na(x$beta.se[i])) && (!is.nan(x$beta.se[i]))) {
        z.value[i] <- x$beta[i]/x$beta.se[i]
        p.value[i] <- 2*pnorm(-abs(z.value[i]),0,1)
      } else {
        z.value[i] <- NA
        p.value[i] <- NA
      } 
    }
    
    aux1 <- nchar(sprintf("%.4f",c(x$lambda,x$xi,x$beta)))
    aux2 <- nchar(sprintf("%.4f",c(x$lambda.se,x$xi.se,x$beta.se)))
    aux3 <- nchar(sprintf("%.4f",z.value))
    p.value2 <- rep(NA,length(x$beta))
    for (i in 1:length(x$beta)) {
      if ((p.value[i] < 0.0001) && (!is.na(p.value[i]))) {
        p.value2[i] <- "< 0.0001"
      } else if (!is.na(p.value[i])) {
        p.value2[i] <- sprintf("%.4f",p.value[i])
      }
    }
    aux4 <- nchar(p.value2)
    text <- c("Estimates","Std. Error","z value","Pr(>|z|)")
    auxc <- nchar(text)
    
    auxc1 <- c(ifelse(max(aux1) > auxc[1],abs(max(aux1)-auxc[1]),0),
               ifelse(max(aux2) > auxc[2],abs(max(aux2)-auxc[2]),0),
               ifelse(max(aux3) > auxc[3],abs(max(aux3)-auxc[3]),0),
               ifelse(max(aux4) > auxc[4],abs(max(aux4)-auxc[4]),0))
    
    cat("         ",rep(" ",auxc1[1]+1),text[1],
        rep(" ",auxc1[2]+1),text[2],
        rep(" ",auxc1[3]+1),text[3],
        rep(" ",auxc1[4]+1),text[4],"\n",sep="")
    auxc <- auxc+auxc1
    cat("  lambda: ",
        rep(ifelse(auxc[1] > aux1[1]," ",""),abs(auxc[1]-aux1[1])),sprintf("%.4f",x$lambda)," ",
        rep(ifelse(auxc[2] > aux2[1]," ",""),abs(auxc[2]-aux2[1])),sprintf("%.4f",x$lambda.se),
        "\n",sep="")
    cat("      xi: ",
        rep(ifelse(auxc[1] > aux1[2]," ",""),abs(auxc[1]-aux1[2])),sprintf("%.4f",x$xi)," ",
        rep(ifelse(auxc[2] > aux2[2]," ",""),abs(auxc[2]-aux2[2])),sprintf("%.4f",x$xi.se),
        "\n",sep="")
    
    if (length(x$beta) == 1) {
      i <- 1
      cat("    beta: ",
          rep(ifelse(auxc[1] > aux1[(i+2)]," ",""),abs(auxc[1]-aux1[(i+2)])),sprintf("%.4f",x$beta[i])," ",
          rep(ifelse(auxc[2] > aux2[(i+2)]," ",""),abs(auxc[2]-aux2[(i+2)])),sprintf("%.4f",x$beta.se[i])," ",
          rep(ifelse(auxc[3] > aux3[i]," ",""),abs(auxc[3]-aux3[i])),sprintf("%.4f",z.value[i])," ",
          rep(ifelse(auxc[4] > aux4[i]," ",""),abs(auxc[4]-aux4[i])),p.value2[i],
          ifelse(!is.na(p.value[i]),ifelse(p.value[i] < 0.0001," ***",
                                    ifelse(p.value[i] < 0.01," **",
                                    ifelse(p.value[i] < 0.05," *",
                                    ifelse(p.value[i] < 0.1," .","")))),""),
          "\n",sep="")
    } else {
      for (i in 1:length(x$beta)) {
        if (i-1 < 10) {
          cat("   beta",i-1,": ",
              rep(ifelse(auxc[1] > aux1[(i+2)]," ",""),abs(auxc[1]-aux1[(i+2)])),sprintf("%.4f",x$beta[i])," ",
              rep(ifelse(auxc[2] > aux2[(i+2)]," ",""),abs(auxc[2]-aux2[(i+2)])),sprintf("%.4f",x$beta.se[i])," ",
              rep(ifelse(auxc[3] > aux3[i]," ",""),abs(auxc[3]-aux3[i])),sprintf("%.4f",z.value[i])," ",
              rep(ifelse(auxc[4] > aux4[i]," ",""),abs(auxc[4]-aux4[i])),p.value2[i],
              ifelse(!is.na(p.value[i]),ifelse(p.value[i] < 0.0001," ***",
                                        ifelse(p.value[i] < 0.01," **",
                                        ifelse(p.value[i] < 0.05," *",
                                        ifelse(p.value[i] < 0.1," .","")))),""),
              "\n",sep="")
        } else {
          cat("  beta",i-1,": ",
              rep(ifelse(auxc[1] > aux1[(i+2)]," ",""),abs(auxc[1]-aux1[(i+2)])),sprintf("%.4f",x$beta[i])," ",
              rep(ifelse(auxc[2] > aux2[(i+2)]," ",""),abs(auxc[2]-aux2[(i+2)])),sprintf("%.4f",x$beta.se[i])," ",
              rep(ifelse(auxc[3] > aux3[i]," ",""),abs(auxc[3]-aux3[i])),sprintf("%.4f",z.value[i])," ",
              rep(ifelse(auxc[4] > aux4[i]," ",""),abs(auxc[4]-aux4[i])),p.value2[i],
              ifelse(!is.na(p.value[i]),ifelse(p.value[i] < 0.0001," ***",
                                        ifelse(p.value[i] < 0.01," **",
                                        ifelse(p.value[i] < 0.05," *",
                                        ifelse(p.value[i] < 0.1," .","")))),""),
              "\n",sep="")
        }
      }
    }
    cat("\n",sep="")
    cat("Summary statistic of the error for the TBS model\n",sep="")
    print(summary(x$error))
    cat("\n",sep="")
    if (length(x$beta) == 1) {
      cat("Estimated quantiles of time event\n",sep="")
      aux1 <- qtbs(c(0.05,0.25,0.5,0.75,0.95),
                   lambda=x$lambda,
                   xi=x$xi,
                   beta=x$beta,
                   dist=x$error.dist)
      aux2 <- nchar(trunc(aux1))
      cat("   ",rep(" ",aux2[1]),"5%   ",rep(" ",aux2[2]),"25%   ",rep(" ",aux2[3]),
          "50%   ",rep(" ",aux2[4]),"75%   ",rep(" ",aux2[5]),"95%\n",sep="")
      cat(sprintf("%.4f",aux1[1])," ",sprintf("%.4f",aux1[2])," ",sprintf("%.4f",aux1[3])," ",
          sprintf("%.4f",aux1[4])," ",sprintf("%.4f",aux1[5]),"\n",sep="")
    }  
    cat("--------------------------------------------------------\n",sep="")
    cat("\n",sep="")
  } else {
    cat("There was not obtained the convergence.\n",sep="")
  }
}

plot.tbs.survreg.mle.best <- function(x, plot.type=NULL, type="l", xlim=NULL, 
                                      ylim=NULL, main=NULL, xlab=NULL, ylab=NULL,
                                      lty=NULL, lwd=NULL, col=NULL, ...) {
  plot.tbs.survreg.mle(x[[x$best.n]], plot.type=plot.type, type=type, xlim=xlim,
                       ylim=ylim, main=main, xlab=xlab, ylab=ylab, lty=lty, lwd=lwd,
                       col=col, ...)
}

plot.tbs.survreg.mle <- function(x, plot.type=NULL, type="l", xlim=NULL, ylim = NULL,
                                 main = NULL, xlab = NULL, ylab = NULL,
                                 lty = NULL, lwd = NULL, col = NULL, ...) {
  if(is.null(plot.type)) {
    if(is.null(x$x)) plot.type <- 'error'
    else plot.type <- 'surv'
  }
  if ((is.null(x$x)) && (! (plot.type %in% c('error'))))
    stop('Invalid plot type for tbs.survreg.mle. It is only possible to use option \'error\' for the estimated model.')
  if(! (plot.type %in% c('surv','hazard','error')))
    stop('Invalid plot type for tbs.survreg.mle. Options are surv, hazard, error.')

  if (x$convergence) {
    h <- 1000
    if (is.null(xlim)) {
      LB <- 0.0001
      UB <- max(x$time)*1.1
      xlim <- c(LB,UB)
    } else {
      LB <- ifelse(xlim[1] <= 0,0.0001,xlim[1])
      UB <- xlim[2]
    }
    axis.t <- seq(LB,UB,(UB-LB)/(h-1))

    k=0
    if(plot.type=='surv') {
      k = 1
      if (is.null(xlab))
        xlab <- "time"
      if (is.null(ylab))
        ylab <- "S(t)"
      if (is.null(main)) {
        main <- "Survival function (MLE)"
      }
    }
    if(plot.type=='hazard') {
      k = 2
      if (is.null(xlab))
        xlab <- "time"
      if (is.null(ylab))
        ylab <- "h(t)"
      if (is.null(main))
        main <- "Hazard function (MLE)"
    }

    if(k > 0) {
      if (attr(x$x,"plot") == 1) {
        if (is.null(lty))
          lty <- 1
        if (is.null(lwd))
          lwd <- 1
        if (is.null(col))
          col <- 1

        if (k == 1) { ### Survival plot
          if (is.null(ylim)) {
            ylim=c(0,1)
          }
          axis.y <- 1-ptbs(axis.t,lambda=x$lambda,xi=x$xi,beta=x$beta,dist=x$error.dist)
          plot(axis.t,axis.y,type=type,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main,
               lwd=lwd,lty=lty,col=col,...)
        } else { ### Hazard plot
          axis.y <- htbs(axis.t,lambda=x$lambda,xi=x$xi,beta=x$beta,dist=x$error.dist)
          if (is.null(ylim)) {
            ylim=c(0,1.1*max(axis.y,na.rm=TRUE))
          }
          plot(axis.t,axis.y,type=type,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main,
               lwd=lwd,lty=lty,col=col,...)
        }
      } else if ((attr(x$x,"plot") == 2) || (attr(x$x,"plot") == 3)) {
        if (is.null(lty)) {
          lty <- seq(1,length(x$x),1)
        } else {
          if (length(lty) != length(x$x)) {
            stop(paste("The 'lty' length must be equal to ",length(x$x),sep=""))
          }
        }
        if (is.null(col)) {
          col <- rep(1,length(x$x))
        } else {
          if (length(col) != length(x$x)) {
            stop(paste("The 'col' length must be equal to ",length(x$x),sep=""))
          }
        }
        if (is.null(lwd)) {
          lwd <- rep(1,length(x$x))
        } else {
          if (length(lwd) != length(x$x)) {
            stop(paste("The 'lwd' length must be equal to ",length(x$x),sep=""))
          }
        }
        axis.y <- matrix(NA,length(axis.t),length(x$x))
        for (i in 1:length(x$x)) {
          if (attr(x$x,"plot") == 2) {
            if (k == 1) { ### Survival plot
              axis.y[,i] <- 1-ptbs(axis.t,lambda=x$lambda,xi=x$xi,beta=x$beta*x$x[i],dist=x$error.dist)
            } else { ### Hazard plot
              axis.y[,i] <- htbs(axis.t,lambda=x$lambda,xi=x$xi,beta=x$beta*x$x[i],dist=x$error.dist)
            }
          } else {
            if (k == 1) { ### Survival plot
              axis.y[,i] <- 1-ptbs(axis.t,lambda=x$lambda,xi=x$xi,
                                   beta=(x$beta[1]+x$beta[2]*x$x[i]),dist=x$error.dist)
            } else { ### Hazard plot
              axis.y[,i] <- htbs(axis.t,lambda=x$lambda,xi=x$xi,
                                 beta=(x$beta[1]+x$beta[2]*x$x[i]),dist=x$error.dist)
            }
          }
        }
        for (i in 1:length(x$x)) {
          if (i == 1) {
            if (k == 1) { ### Survival plot
              if (is.null(ylim)) {
                ylim=c(0,1)
              }
              plot(axis.t,axis.y[,i],type=type,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main,
                   lty=lty[i],lwd=lwd[i],col=col[i], ...)
            } else {
              if (is.null(ylim)) {
                ylim=c(0,1.1*max(axis.y,na.rm=TRUE))
              }
              plot(axis.t,axis.y[,i],type=type,xlab=xlab,ylab=ylab,main=main,
                   ylim=ylim,xlim=xlim,lty=lty[i],lwd=lwd[i],col=col[i], ...)
            }
          } else {
            lines(axis.t,axis.y[,i],lty=lty[i],lwd=lwd[i],col=col[i])
          }
        }
      }
    }
    if(plot.type=='error') {
      ## Histogram for error distribution:
      if (is.null(xlab))
        xlab <- "error"
      if (is.null(ylab))
        ylab <- "Density"
      if (is.null(main))
        main <- "Error Distribution (MLE)"
      bound <- max(abs(c(min(x$error),max(x$error))))
      hist(x$error,probability=TRUE,xlim=c(-bound,bound),main=main,xlab=xlab,ylab=ylab)
      lines(seq(-bound,bound,2*bound/(h-1)),
            x$error.dist$d(seq(-bound,bound,2*bound/(h-1)),xi=x$xi))
    }
    if(plot.type=='qqplot') {
      ## Q-Q plot for error distribution:
      qqplot(x$error.dist$q(ppoints(length(x$error)), xi=x$xi),x$error,
             main = expression("Q-Q plot for error"),xlab="Theoretical Quantiles",ylab="Sample Quantiles")
      qqline(x$error, distribution = function(p) x$error.dist$q(p, xi=x$xi),
             probs = c(0.25, 0.75), col = 2, lwd=2)
    }
  } else {
    cat("Convergence has not been obtained for this tbs.survreg.mle.\n",sep="")
  }
}


lines.tbs.survreg.mle.best <- function(x, plot.type='surv', lty=NULL, lwd=NULL, col=NULL, ...) {
  lines.tbs.survreg.mle(x[[x$best.n]], plot.type=plot.type, lty=lty, lwd=lwd, col=col, ...)
}

lines.tbs.survreg.mle <- function(x, plot.type="surv", lty=NULL, lwd=NULL, col=NULL, ...) {
  if(! (plot.type %in% c('surv','hazard')))
    stop('Invalid plot type for tbs.survreg.mle. Options are surv, hazard.')
  if (is.null(x$x))
    stop('Invalid plot type for tbs.survreg.mle. It is not possible to draw a plot for the estimated model.')

  if (x$convergence) {
    h <- 1000

    LB <- 0.0001
    UB <- max(x$time)*1.1
    axis.t <- seq(LB,UB,(UB-LB)/(h-1))

    if (attr(x$x,"plot") == 1) {
      if (is.null(lty))
        lty <- 1
      if (is.null(lwd))
        lwd <- 1
      if (is.null(col))
        col <- 1

      if (plot.type=='surv') { ### Survival plot
        axis.y <- 1-ptbs(axis.t,lambda=x$lambda,xi=x$xi,beta=x$beta,dist=x$error.dist)
      } else { ### Hazard plot
        axis.y <- htbs(axis.t,lambda=x$lambda,xi=x$xi,beta=x$beta,dist=x$error.dist)
      }
      lines(axis.t,axis.y, lty=lty, lwd=lwd, col=col, ...)
    } else if ((attr(x$x,"plot") == 2) || (attr(x$x,"plot") == 3)) {
      if (is.null(lty)) {
        lty <- seq(1,length(x$x),1)
      } else {
        if (length(lty) != length(x$x)) {
          stop(paste("The 'lty' length must be equal to ",length(x$x),sep=""))
        }
      }
      if (is.null(col)) {
        col <- rep(1,length(x$x))
      } else {
        if (length(col) != length(x$x)) {
          stop(paste("The 'col' length must be equal to ",length(x$x),sep=""))
        }
      }
      if (is.null(lwd)) {
        lwd <- rep(1,length(x$x))
      } else {
        if (length(lwd) != length(x$x)) {
          stop(paste("The 'lwd' length must be equal to ",length(x$x),sep=""))
        }
      }
      axis.y <- matrix(NA,length(axis.t),length(x$x))
      for (i in 1:length(x$x)) {
        if (attr(x$x,"plot") == 2) {
          if (plot.type=='surv') { ### Survival plot
            axis.y[,i] <- 1-ptbs(axis.t,lambda=x$lambda,xi=x$xi,beta=x$beta*x$x[i],dist=x$error.dist)
          } else { ### Hazard plot
            axis.y[,i] <- htbs(axis.t,lambda=x$lambda,xi=x$xi,beta=x$beta*x$x[i],dist=x$error.dist)
          }
        } else {
          if (plot.type=='surv') { ### Survival plot
            axis.y[,i] <- 1-ptbs(axis.t,lambda=x$lambda,xi=x$xi,
                                 beta=(x$beta[1]+x$beta[2]*x$x[i]),dist=x$error.dist)
          } else { ### Hazard plot
            axis.y[,i] <- htbs(axis.t,lambda=x$lambda,xi=x$xi,
                               beta=(x$beta[1]+x$beta[2]*x$x[i]),dist=x$error.dist)
          }
        }
      }
      for (i in 1:length(x$x)) {
        lines(axis.t,axis.y[,i],lty=lty[i],lwd=lwd[i],col=col[i], ...)
      }
    }
  } else {
    cat("Convergence has not been obtained for this tbs.survreg.mle.\n",sep="")
  }
}

