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

## TBS estimation using a Bayesian approach. The lack of a closed form
## solution forces us to tackle the problem by MCMC.
## max.time is a time limit in minutes (<= 0 means no limit)
## formula is specification containing a Surv model with right-censored data as in the package survival.
## dist defines the error distribution by name ("norm", "doubexp", "t", "cauchy", "logistic") or using
## a call from dist.error, or even created directly by the user (see dist.error.r for details)
## guess.beta, guess.lambda, guess.xi are initial value of the Markov Chain (beta has to have the same
## of elements as covariates, lambda and xi are scalars.
## burn-in is the number of firsts samples of posterior to not use, and jump is the number of jump
## between each sample of posterior to avoid the problem of auto-correlation between the samples.
## size is the final sample size of the posterior. scale is a parameter of the `metrop' function which
## controls the acceptance rate. prior.mean and prior.sd define the parameters of the normal prior for
## the MCMC (by default, they are equal to 5 and 5).
## accept: fraction of Metropolis proposals accepted.
tbs.survreg.be <- function(formula,dist=dist.error("norm"),max.time=-1,
                           guess.beta=NULL,guess.lambda=1,guess.xi=1,
                           burn=1000,jump=2,size=500,scale=0.1,
                           prior.mean=NULL,prior.sd=NULL,
                           seed=1234) {
  if(is.character(dist)) dist=dist.error(dist)
  initial.time <- .gettime()
  set.seed(seed)
  if(max.time <= 0) {
    ## user didn't define a timeout, so we set to a large number
    max.time <- 1e10
  }
  ## check the class of formula
  if (attributes(formula)$class != "formula")
    stop("A formula argument is required")

  ## record the call arguments
  Call  <- match.call()
  ## read the information from within the formula to populate the required variables
  mf <- model.frame(formula=formula)
  x <- model.matrix(attr(mf, "terms"), data=mf)
  y <- model.response(mf)
  time <- y[,1]
  delta <- y[,2]
  x.k   <- dim(x)[2]
  n     <- dim(x)[1]
  ## check if delta is an indicator function
  if (any((delta != 0) & (delta != 1)))  {
    stop("It is only accepted uncensored or right censored data")
  }
  if (is.null(guess.beta)) {
    guess.beta <- rep(0,x.k)
  }

  out <- NULL
  out$call <- Call
  if (is.matrix(x)) {
#    out$x <- x[order(time),]
    x <- x[order(time),]
  } else {
#    out$x <- x[order(time)]
    x <- x[order(time)]
  }
  out$delta <- delta[order(time)]
  out$time  <- time[order(time)]
#  x     <- out$x
  delta <- out$delta
  time  <- out$time
##  out$time <- time[delta == 1]

  ## perform a series of verifications for the given arguments of the function
  if (length(guess.lambda) != 1)
    stop("guess.lambda is not a scalar")
  if (guess.lambda <= 0)
    stop("guess.lambda must be a positive number")
  if (length(guess.xi) != 1)
    stop("guess.xi is not a scalar")
  if (guess.xi <= 0)
    stop("guess.xi must be a positive number")
  if (length(guess.beta) != x.k)
    stop("guess.beta length is not in accordance with the model specification")
  guess <- c(guess.lambda,guess.xi,guess.beta)

  if ((!is.integer(burn)) && (burn <= 0)) 
    stop("burn must be a positive integer number")
  if ((!is.integer(jump)) && (jump <= 0)) 
    stop("jump must be a positive integer number")
  if ((!is.integer(size)) && (size <= 0)) 
    stop("size must be a positive integer number")
  if (scale <= 0)
    stop("scale must be a positive number")

  if (!is.null(x)) {
    if (is.matrix(x)) {
      if (length(time) != length(x[,1]))
        stop("length of time is different of length of x")
    } else {
      if (length(time) != length(x))
        stop("length of time is different of length of x")
      x <- matrix(x,length(x),1)
    }
    if(length(beta) > 1)
      beta <- matrix(beta,length(beta),1)
  } else {
    x <- matrix(1,length(time),1)
  }

  if (is.null(prior.mean)) {
    prior.mean <- 0
  } else {
    if (!is.vector(prior.mean))
      stop("mean is not a vector/scalar")
    if ((length(prior.mean) != 1) || (length(prior.mean) != length(guess[3:length(guess)]))) {
      stop(paste("length mean is different of 1 or ",length(guess[3:length(guess)]),sep=""))
    }
  }
  if (is.null(prior.sd)) {
    prior.sd <- 100
  } else {
    if (!is.vector(prior.sd))
      stop("sd is not a vector/scalar")
    if ((length(prior.sd) != 1) || (length(prior.sd) != length(guess[3:length(guess)]))) {
      stop(paste("length sd is different of 1 or ",length(par[3:length(guess)]),sep=""))
    }
    if (prior.sd <= 0)
      stop("prior.sd must be a positive number")
  }

  ## call the Metropolis algorithm for MCMC
  chain <- try(withTimeout(metrop(obj=.logpost,initial=guess,time=time,delta=delta,dist=dist,x=x,
                                      mean=prior.mean,sd=prior.sd,
                                      nbatch=(size-1)*jump+burn,blen=1,nspac=1,scale=scale),
                           timeout=max.time*60,onTimeout="error"),silent=TRUE)
  if(length(class(chain)) == 1 && class(chain) == "try-error") {
    stop("Time limit exceeded")
  }
  out$post <- chain$batch[seq(burn,length(chain$batch[,1]),jump),]
  out$accept <- chain$accept

  # evaluating the point estimates
  if (x.k != 1) {
    par    <- c(mean(out$post[,1]),mean(out$post[,2]),apply(out$post[,3:length(out$post[1,])],2,mean))
    par.sd <- c(sd(out$post[,1]),    sd(out$post[,2]),apply(out$post[,3:length(out$post[1,])],2,sd))
  } else { 
    par    <- c(mean(out$post[,1]),mean(out$post[,2]),mean(out$post[,3:length(out$post[1,])]))
    par.sd <- c(sd(out$post[,1]),    sd(out$post[,2]),  sd(out$post[,3:length(out$post[1,])]))
  }
  out$lambda <- par[1]
  out$xi     <- par[2]
  out$beta   <- par[3:length(par)]
  out$lambda.sd <- par.sd[1]
  out$xi.sd     <- par.sd[2]
  out$beta.sd   <- par.sd[3:length(par.sd)]

  text <- c("lambda","xi")
  if (length(out$beta) == 1) {
    text <- c(text,"beta")
  } else {
    for (i in 1:length(out$beta))
      text <- c(text,paste("beta",i-1,sep=""))
  }
  colnames(out$post) <- text
  # evaluating the interval estimates
  par.HPD   <- cbind(c(HPDinterval(as.mcmc(out$post[,1]),0.95)),
                     c(HPDinterval(as.mcmc(out$post[,2]),0.95)))
  for (i in 3:length(out$post[1,])) {
    par.HPD <- cbind(par.HPD,
                     c(HPDinterval(as.mcmc(out$post[,i]),0.95)))
  }
  out$lambda.HPD <- par.HPD[,1]
  out$xi.HPD     <- par.HPD[,2]
  out$beta.HPD   <- par.HPD[,3:length(par.HPD[1,])]

  # evaluating DIC
  aux.loglik  <- rep(0,size)
  for (j in 1:size)
    aux.loglik[j] <- c(.lik.tbs(out$post[j,],time,delta,dist,x))
  loglik <- mean(-2*aux.loglik)
  rm(aux.loglik)
  out$DIC <- 2*loglik+2*.lik.tbs(par,time,delta,dist,x)

  # evaluating error of the model
  aux.error    <- matrix(0,length(time),size)
  for (j in 1:size) {
      if (x.k != 1) {
        aux.error[,j]    <- (.g.lambda(log(time),out$post[j,1])-
           .g.lambda(x%*%matrix(out$post[j,3:length(out$post[1,])],length(out$post[j,3:length(out$post[1,])]),1),out$post[j,1]))
      } else {
        aux.error[,j]    <- (.g.lambda(log(time),out$post[j,1])-
                             .g.lambda(out$post[j,3]*x,out$post[j,1]))
      }
  }
  error <- rep(0,length(time))
  for (i in 1:length(time))
    error[i] <- mean(aux.error[i,])
  rm(aux.error)
  out$error <- error
##  out$error <- error[delta == 1]
  out$error.dist <- dist

  aux <- .test.tbs(out$lambda,out$xi,out$beta,x,time,type="d")
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
  # time spent to do BE
  out$run.time <- .gettime() - initial.time

  class(out) <- "tbs.survreg.be"
  return(out)
}

print.tbs.survreg.be <- function(x, ...) {
#  print(x$call)
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
  cat(" TBS model with ",text.dist," error distribution (BE).\n",sep="")
  cat("\n",sep="")

  aux1 <- nchar(sprintf("%.4f",c(x$lambda,x$xi,x$beta)))
  aux2 <- nchar(sprintf("%.4f",c(x$lambda.sd,x$xi.sd,x$beta.sd)))
  text <- c("Estimates","Std. Error")
  auxc <- nchar(text)
  auxc1 <- c(ifelse(max(aux1) > auxc[1],abs(max(aux1)-auxc[1]),0),
             ifelse(max(aux2) > auxc[2],abs(max(aux2)-auxc[2]),0))
  auxc <- auxc+auxc1

  cat("         ",rep(" ",auxc1[1]+1),text[1],
                  rep(" ",auxc1[2]+1),text[2],"\n",sep="")
  cat("  lambda: ",
      rep(ifelse(auxc[1] > aux1[1]," ",""),abs(auxc[1]-aux1[1])),sprintf("%.4f",x$lambda)," ",
      rep(ifelse(auxc[2] > aux2[1]," ",""),abs(auxc[2]-aux2[1])),sprintf("%.4f",x$lambda.sd),
      "\n",sep="")
  cat("      xi: ",
      rep(ifelse(auxc[1] > aux1[2]," ",""),abs(auxc[1]-aux1[2])),sprintf("%.4f",x$xi)," ",
      rep(ifelse(auxc[2] > aux2[2]," ",""),abs(auxc[2]-aux2[2])),sprintf("%.4f",x$xi.sd),
      "\n",sep="")
  if (length(x$beta) == 1) {
    cat("    beta: ",
        rep(ifelse(auxc[1] > aux1[3]," ",""),abs(auxc[1]-aux1[3])),sprintf("%.4f",x$beta)," ",
        rep(ifelse(auxc[2] > aux2[3]," ",""),abs(auxc[2]-aux2[3])),sprintf("%.4f",x$beta.sd)," ",
        "\n",sep="")
  } else {
    for (i in 1:length(x$beta)) {
      if (i-1 < 10) {
        cat("   beta",i-1,": ",
            rep(ifelse(auxc[1] > aux1[(i+2)]," ",""),abs(auxc[1]-aux1[(i+2)])),sprintf("%.4f",x$beta[i])," ",
            rep(ifelse(auxc[2] > aux2[(i+2)]," ",""),abs(auxc[2]-aux2[(i+2)])),sprintf("%.4f",x$beta.sd[i])," ",
            "\n",sep="")
      } else {
        cat("  beta",i-1,": ",
            rep(ifelse(auxc[1] > aux1[(i+2)]," ",""),abs(auxc[1]-aux1[(i+2)])),sprintf("%.4f",x$beta[i])," ",
            rep(ifelse(auxc[2] > aux2[(i+2)]," ",""),abs(auxc[2]-aux2[(i+2)])),sprintf("%.4f",x$beta.sd[i])," ",
            "\n",sep="")
      }
    }
  }
  cat("\n",sep="")
  cat("     DIC: ",sprintf("%.4f",x$DIC),"\n",sep="")
  cat("\n",sep="")
  cat("Run time: ", sprintf("%.2f",x$run.time)," min \n",sep="")
  cat("--------------------------------------------------------\n",sep="")
  cat("\n",sep="")
}

summary.tbs.survreg.be <- function(object, ...) {
  x=object
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
  cat(" TBS model with ",text.dist," error distribution (BE).\n",sep="")
  cat("\n",sep="")

  auxt1 <- sprintf("%.4f",c(x$lambda,x$xi,x$beta))
  aux1  <- nchar(auxt1)
  auxt2 <- sprintf("%.4f",c(x$lambda.sd,x$xi.sd,x$beta.sd))
  aux2  <- nchar(auxt2)
  if (length(x$beta) ==  1) {
    auxt3 <- sprintf("%.4f",c(x$lambda.HPD[1],x$xi.HPD[1],x$beta.HPD[1]))
    aux3  <- nchar(auxt3)
    auxt4 <- sprintf("%.4f",c(x$lambda.HPD[2],x$xi.HPD[2],x$beta.HPD[2]))
    aux4  <- nchar(auxt4)
  } else {
    auxt3 <- sprintf("%.4f",c(x$lambda.HPD[1],x$xi.HPD[1],x$beta.HPD[1,]))
    aux3  <- nchar(auxt3)
    auxt4 <- sprintf("%.4f",c(x$lambda.HPD[2],x$xi.HPD[2],x$beta.HPD[2,]))
    aux4  <- nchar(auxt4)
  }

  text <- c("     Mean"," Std. Dev.","95% HPD CI")
  text2 <- c("  lambda: ","      xi: ")
  if (length(x$beta) == 1) {
    text2 <- c(text2,"    beta: ")
  } else {
    for (i in 1:length(x$beta)) {
      if (i-1 < 10) {
        text2 <- c(text2,paste("   beta",i-1,": ",sep=""))
      } else {
        text2 <- c(text2,paste("  beta",i-1,": ",sep=""))
      }
    }
  }
  auxc <- nchar(text)

  auxc2    <- max(aux1,auxc[1]) -auxc[1]
  auxc2[2] <- max(aux2,auxc[2]) -auxc[2]
  auxc2[3] <- max(aux3+aux4+3,auxc[1]) -auxc[3]
  auxc2 <- ifelse(auxc2 < 0,0,auxc2)
  auxc[3] <- max(aux3)
  auxc[4] <- max(aux4)

  cat("         ",rep(" ",auxc2[1]+1),text[1],
                  rep(" ",auxc2[2]+1),text[2],
                  rep(" ",auxc2[3]+1),text[3],"\n",sep="")
  for (i in 1:length(text2)) {
    cat(text2[i],
        rep(ifelse(auxc[1] > aux1[i]," ",""),abs(auxc[1]-aux1[i])),auxt1[i]," ",
        rep(ifelse(auxc[2] > aux2[i]," ",""),abs(auxc[2]-aux2[i])),auxt2[i]," (",
        rep(ifelse(auxc[3] > aux3[i]," ",""),abs(auxc[3]-aux3[i])),auxt3[i],", ",
        rep(ifelse(auxc[4] > aux4[i]," ",""),abs(auxc[4]-aux4[i])),auxt4[i],")",
        "\n",sep="")
  }
  cat("\n",sep="")
  cat("Summary statistic of the posterior mean of the error\n")
  cat("for the TBS model:\n",sep="")
  print(summary(x$error))
  cat("\n",sep="")
  if (length(x$beta) == 1) {
    cat("Estimated quantiles of time event:\n",sep="")
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
}

plot.tbs.survreg.be <- function(x, plot.type='surv', HPD=TRUE, HPD.alpha=0.95, type="l",
                                xlim=NULL, ylim = NULL, main = NULL, xlab = NULL,
                                ylab = NULL, lty = NULL, lwd = NULL, col = NULL,
                                lty.HPD = NULL, lwd.HPD = NULL, col.HPD = NULL, ...) {
  if(! (plot.type %in% c('surv','hazard','error','auto','ts')))
    stop('Invalid plot type for tbs.survreg.be. Options are surv, hazard, error, auto, ts.')
  if ((is.null(x$x)) && (! (plot.type %in% c('error','auto','ts'))))
    stop('Invalid plot type for tbs.survreg.mle. It is only possible to use the options \'error\',\'auto\' or \'ts\' for the estimated model.')
  if ((HPD.alpha >= 1) || (HPD.alpha <= 0))
    stop('HPD.alpha must be between 0 and 1.')

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
      main <- "Survival function (BE)"
    }
  }
  if(plot.type=='hazard') {
    k = 2
    if (is.null(xlab))
      xlab <- "time"
    if (is.null(ylab))
      ylab <- "h(t)"
    if (is.null(main))
      main <- "Hazard function (BE)"
  }

  if(k > 0) {
    if (attr(x$x,"plot") == 1) {
      if (is.null(lty))
        lty <- 1
      if (is.null(lwd))
        lwd <- 1
      if (is.null(col))
        col <- 1

      axis.y <- matrix(NA,length(axis.t),length(x$post[,1]))
      if (k == 1) { ### Survival plot
        if (is.null(ylim)) {
          ylim=c(0,1)
        }

        for (j in 1:length(x$post[,1])) {
          axis.y[,j] <- 1-ptbs(axis.t,lambda=x$post[j,1],xi=x$post[j,2],
                               beta=x$post[j,3],dist=x$error.dist)
        }
        plot(axis.t,apply(axis.y,1,mean),type=type,xlim=xlim,ylim=ylim,xlab=xlab,
             ylab=ylab,main=main,lwd=lwd,lty=lty,col=col,...)
      } else { ### Hazard plot
        for (j in 1:length(x$post[,1])) {
          axis.y[,j] <- htbs(axis.t,lambda=x$post[j,1],xi=x$post[j,2],
                             beta=x$post[j,3],dist=x$error.dist)
        }
        if (is.null(ylim)) {
          ylim=c(0,1.1*max(apply(axis.y,1,quantile,probs=(1-(1-HPD.alpha)/2)),na.rm=TRUE))
        }
        plot(axis.t,apply(axis.y,1,mean),type=type,xlim=xlim,ylim=ylim,xlab=xlab,
             ylab=ylab,main=main,lwd=lwd,lty=lty,col=col,...)
      }
      if (HPD) {
        if (is.null(lty.HPD))
          lty.HPD <- 1
        if (is.null(lwd.HPD))
          lwd.HPD <- 1
        if (is.null(col.HPD))
          col.HPD <- "gray50"
        lines(axis.t,apply(axis.y,1,quantile,probs=(1-HPD.alpha)/2),col=col.HPD,lty=lty.HPD,lwd=lwd.HPD, ...)
        lines(axis.t,apply(axis.y,1,quantile,probs=1-(1-HPD.alpha)/2),col=col.HPD,lty=lty.HPD,lwd=lwd.HPD, ...)
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

      axis.y <- array(NA,c(length(axis.t),length(x$post[,1]),length(x$x)))
      for (i in 1:length(x$x)) {
        if (attr(x$x,"plot") == 2) {
          if (k == 1) { ### Survival plot
            for (j in 1:length(x$post[,1])) {
              axis.y[,j,i] <- 1-ptbs(axis.t,lambda=x$post[j,1],xi=x$post[j,2],
                                   beta=x$post[j,3]*x$x[i],dist=x$error.dist)
            }
          } else { ### Hazard plot
            for (j in 1:length(x$post[,1])) {
              axis.y[,j,i] <- htbs(axis.t,lambda=x$post[j,1],xi=x$post[j,2],
                                 beta=x$post[j,3]*x$x[i],dist=x$error.dist)
            }
          }
        } else {
          if (k == 1) { ### Survival plot
            for (j in 1:length(x$post[,1])) {
              axis.y[,j,i] <- 1-ptbs(axis.t,lambda=x$post[j,1],xi=x$post[j,2],
                                     beta=(x$post[j,3]+x$post[j,4]*x$x[i]),dist=x$error.dist)
            }
          } else { ### Hazard plot
            for (j in 1:length(x$post[,1])) {
              axis.y[,j,i] <- htbs(axis.t,lambda=x$post[j,1],xi=x$post[j,2],
                                   beta=(x$post[j,3]+x$post[j,4]*x$x[i]),dist=x$error.dist)
            }
          }
        }
      }
      aux.HPD <- matrix(NA,length(axis.t),length(x$x))
      for (i in 1:length(x$x)) {
        aux.HPD[,i] <- apply(axis.y[,,i],1,quantile,probs=(1-(1-HPD.alpha)/2))
      }
      if (HPD) {
        if (is.null(lty.HPD)) {
          lty.HPD <- seq(1,length(x$x),1)
        } else {
          if (length(lty.HPD) != length(x$x)) {
            stop(paste("The 'lty.HPD' length must be equal to ",length(x$x),sep=""))
          }
        }
        if (is.null(col.HPD)) {
          col.HPD <- rep("gray50",length(x$x))
        } else {
          if (length(col.HPD) != length(x$x)) {
            stop(paste("The 'col.HPD' length must be equal to ",length(x$x),sep=""))
          }
        }
        if (is.null(lwd.HPD)) {
          lwd.HPD <- rep(1,length(x$x))
        } else {
          if (length(lwd.HPD) != length(x$x)) {
            stop(paste("The 'lwd.HPD' length must be equal to ",length(x$x),sep=""))
          }
        }
      }
      for (i in 1:length(x$x)) {
        if (i == 1) {
          if (k == 1) { ### Survival plot
            if (is.null(ylim)) {
              ylim=c(0,1)
            }
            plot(axis.t,apply(axis.y[,,i],1,mean),type=type,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main,
                 lty=lty[i],lwd=lwd[i],col=col[i], ...)
          } else { ### Hazard plot
            if (is.null(ylim)) {
              ylim=c(0,1.1*max(aux.HPD,na.rm=TRUE))
            }
            plot(axis.t,apply(axis.y[,,i],1,mean),type=type,xlab=xlab,ylab=ylab,main=main,
                 ylim=ylim,xlim=xlim,lty=lty[i],lwd=lwd[i],col=col[i], ...)
          }
        } else {
          lines(axis.t,apply(axis.y[,,i],1,mean),lty=lty[i],lwd=lwd[i],col=col[i])
        }
        if (HPD) {
          lines(axis.t,apply(axis.y[,,i],1,quantile,probs=(1-HPD.alpha)/2),col=col.HPD[i],lty=lty.HPD[i],lwd=lwd.HPD[i], ...)
          lines(axis.t,aux.HPD[,i],col=col.HPD[i],lty=lty.HPD[i],lwd=lwd.HPD[i], ...)
        }
      }
    }
  }
  if(plot.type=='error') {
    ## Histogram for error distribution:
    if (!exists("xlab"))
      xlab <- "error"
    if (!exists("ylab"))
      ylab <- "Density"
    if (!exists("main"))
      main <- "Posterior Mean of Error Distribution (BE)"

    bound <- max(abs(c(min(x$error),max(x$error))))
    hist(x$error,probability=TRUE,xlim=c(-bound,bound),main=main,xlab=xlab,ylab=ylab)
    axis.x <- seq(-bound,bound,2*bound/(h-1))
    axis.y <- matrix(NA,length(axis.x),length(x$post[,1]))
    for (j in 1:length(x$post[,1])) {
      axis.y[,j] <- x$error.dist$d(axis.x,xi=x$post[j,2])
    }

    if (is.null(lty))
      lty <- 1
    if (is.null(lwd))
      lwd <- 1
    if (is.null(col))
      col <- 1
    lines(axis.x,apply(axis.y,1,mean),col=col,lty=lty,lwd=lwd, ...)

    if (HPD) {
      if (is.null(lty.HPD))
        lty.HPD <- 1
      if (is.null(lwd.HPD))
        lwd.HPD <- 1
      if (is.null(col.HPD))
        col.HPD <- "gray50"
      lines(axis.x,apply(axis.y,1,quantile,probs=(1-HPD.alpha)/2),col=col.HPD,lty=lty.HPD,lwd=lwd.HPD, ...)
      lines(axis.x,apply(axis.y,1,quantile,probs=(1-(1-HPD.alpha)/2)),col=col.HPD,lty=lty.HPD,lwd=lwd.HPD, ...)
    }
  }
  if (FALSE) {
    ## Q-Q plot for error distribution:
    axis.x <- matrix(NA,length(x$error),length(x$post[,1]))
    for (j in 1:length(x$post[,1])) {
      axis.x[,j] <- qqplot(x$error.dist$q(ppoints(length(x$error)), xi=x$post[j,2]),x$error,plot.it=FALSE)$x
    }
    axis.y <- qqplot(x$error.dist$q(ppoints(length(x$error)), xi=x$xi),x$error,plot.it=FALSE)$y
    plot(apply(axis.x,1,mean),axis.y,
         main = expression("Q-Q plot for error"),xlab="Theoretical Quantiles",ylab="Sample Quantiles")
    qqline(x$error, distribution = function(p) {
                                     fn.aux <- function(p) {
                                       aux <- rep(NA,length(x$post[,1]))
                                       for (j in 1:length(x$post[,1]))
                                         aux[j] <- x$error.dist$q(p, xi=x$post[j,2])
                                       return(mean(aux))
                                     }
                                     aux <- Vectorize(fn.aux,"p")
                                     return(aux(p))
                                   },
           probs = c(0.25, 0.75), col = 2, lwd=2)
  }

  if(plot.type=='auto') {
    ## autocorrelation plot
    par(ask=TRUE)
    acf(x$post[,1],main=expression(lambda))
    acf(x$post[,2],main=expression(xi))
    if (length(x$post[1,]) == 3) {
      acf(x$post[,3],main=expression(beta))
    } else {
      for (i in 3:length(x$post[1,])) {
        acf(x$post[,i],main=bquote(beta[.(i-3)]))
      }
    }
    par(ask=FALSE)
  }

  if(plot.type=='ts') {
    ## time series plot of the MCMC
    par(ask=TRUE)
    plot(ts(x$post[,1]),ylab=expression(lambda),xlab="i")
    plot(ts(x$post[,2]),ylab=expression(xi),xlab="i")
    if (length(x$post[1,]) == 3) {
      plot(ts(x$post[,3]),ylab=expression(beta),xlab="i")
    } else {
      for (i in 3:length(x$post[1,])) {
        plot(ts(x$post[,i]),ylab=bquote(beta[.(i-3)]),xlab="i")
      }
    }
    par(ask=FALSE)
  }
}

lines.tbs.survreg.be <- function(x, plot.type='surv', HPD=TRUE, HPD.alpha = 0.95,
                                 lty = NULL, lwd = NULL, col = NULL,
                                 lty.HPD = NULL, lwd.HPD = NULL, col.HPD = NULL, ...) {
  if(! (plot.type %in% c('surv','hazard')))
    stop('Invalid plot type for tbs.survreg.be. Options are surv and hazard.')
  if (is.null(x$x))
    stop('Invalid plot type for tbs.survreg.mle. It is not possible to draw a plot for the estimated model.')
  if ((HPD.alpha >= 1) || (HPD.alpha <= 0))
    stop('HPD.alpha must be between 0 and 1.')

  h <- 1000

  LB <- 0.0001
  UB <- max(x$time)*1.1
  axis.t <- seq(LB,UB,(UB-LB)/(h-1))

  if (plot.type == "surv") {
    k <- 1 
  } else if (plot.type == "hazard") {
    k <- 2
  }

  if (attr(x$x,"plot") == 1) {
    if (is.null(lty))
      lty <- 1
    if (is.null(lwd))
      lwd <- 1
    if (is.null(col))
      col <- 1

    if (HPD) {
      if (is.null(lty.HPD))
        lty.HPD <- 1
      if (is.null(lwd.HPD))
        lwd.HPD <- 1
      if (is.null(col.HPD))
        col.HPD <- "gray50"
    }

    axis.y <- matrix(NA,length(axis.t),length(x$post[,1]))
    if (k == 1) { ### Survival plot
      for (j in 1:length(x$post[,1])) {
        axis.y[,j] <- 1-ptbs(axis.t,lambda=x$post[j,1],xi=x$post[j,2],
                             beta=x$post[j,3],dist=x$error.dist)
      }
    } else { ### Hazard plot
      for (j in 1:length(x$post[,1])) {
        axis.y[,j] <- htbs(axis.t,lambda=x$post[j,1],xi=x$post[j,2],
                           beta=x$post[j,3],dist=x$error.dist)
      }
    }
    lines(axis.t,apply(axis.y,1,mean),lwd=lwd,lty=lty,col=col,...)
    if (HPD) {
      lines(axis.t,apply(axis.y,1,quantile,probs=(1-HPD.alpha)/2),col=col.HPD,lty=lty.HPD,lwd=lwd.HPD, ...)
      lines(axis.t,apply(axis.y,1,quantile,probs=(1-(1-HPD.alpha)/2)),col=col.HPD,lty=lty.HPD,lwd=lwd.HPD, ...)
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

    if (HPD) {
      if (is.null(lty.HPD)) {
        lty.HPD <- seq(1,length(x$x),1)
      } else {
        if (length(lty.HPD) != length(x$x)) {
          stop(paste("The 'lty.HPD' length must be equal to ",length(x$x),sep=""))
        }
      }
      if (is.null(col.HPD)) {
        col.HPD <- rep("gray50",length(x$x))
      } else {
        if (length(col.HPD) != length(x$x)) {
          stop(paste("The 'col.HPD' length must be equal to ",length(x$x),sep=""))
        }
      }
      if (is.null(lwd.HPD)) {
        lwd.HPD <- rep(1,length(x$x))
      } else {
        if (length(lwd.HPD) != length(x$x)) {
          stop(paste("The 'lwd.HPD' length must be equal to ",length(x$x),sep=""))
        }
      }
    }

    axis.y <- array(NA,c(length(axis.t),length(x$post[,1]),length(x$x)))
    for (i in 1:length(x$x)) {
      if (attr(x$x,"plot") == 2) {
        if (k == 1) { ### Survival plot
          for (j in 1:length(x$post[,1])) {
            axis.y[,j,i] <- 1-ptbs(axis.t,lambda=x$post[j,1],xi=x$post[j,2],
                                 beta=x$post[j,3]*x$x[i],dist=x$error.dist)
          }
        } else { ### Hazard plot
          for (j in 1:length(x$post[,1])) {
            axis.y[,j,i] <- htbs(axis.t,lambda=x$post[j,1],xi=x$post[j,2],
                               beta=x$post[j,3]*x$x[i],dist=x$error.dist)
          }
        }
      } else {
        if (k == 1) { ### Survival plot
          for (j in 1:length(x$post[,1])) {
            axis.y[,j,i] <- 1-ptbs(axis.t,lambda=x$post[j,1],xi=x$post[j,2],
                                   beta=(x$post[j,3]+x$post[j,4]*x$x[i]),dist=x$error.dist)
          }
        } else { ### Hazard plot
          for (j in 1:length(x$post[,1])) {
            axis.y[,j,i] <- htbs(axis.t,lambda=x$post[j,1],xi=x$post[j,2],
                                 beta=(x$post[j,3]+x$post[j,4]*x$x[i]),dist=x$error.dist)
          }
        }
      }
    }
    for (i in 1:length(x$x)) {
      lines(axis.t,apply(axis.y[,,i],1,mean),lty=lty[i],lwd=lwd[i],col=col[i])
      if (HPD) {
        lines(axis.t,apply(axis.y[,,i],1,quantile,probs=(1-HPD.alpha)/2),col=col.HPD[i],lty=lty.HPD[i],lwd=lwd.HPD[i], ...)
        lines(axis.t,apply(axis.y[,,i],1,quantile,probs=(1-(1-HPD.alpha)/2)),col=col.HPD[i],lty=lty.HPD[i],lwd=lwd.HPD[i], ...)
      }
    }
  }
}


