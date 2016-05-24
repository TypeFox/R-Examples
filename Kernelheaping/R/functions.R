#' @importFrom grDevices rainbow
#' @importFrom graphics box contour image legend lines par plot
#' @importFrom stats aggregate density dnorm optim pnorm qnorm quantile rgamma runif rnorm
#' @importFrom sp point.in.polygon
#' @importFrom plyr round_any count
Rprx=function(rounds, Bias, RR, beta, gridx, ranefvalues, roundvalues){
  rprobs=rprobsRR2(RR,beta,gridx,ranefvalues)
  postMatrix=matrix(0,ncol=ncol(rprobs),nrow=2*nrow(rprobs))
  for(i in 1:length(rounds)){
    modulo=gridx%%rounds[i]
    lower=which(modulo<0.5*rounds[i])
    higher=which(modulo>0.5*rounds[i])
    postMatrix[i,lower]=Bias*rprobs[i,lower]
    postMatrix[i+length(rounds),higher]=(1-Bias)*rprobs[i,higher]
  }
  if(min(roundvalues)==0){postMatrix[c(1,1+length(rounds)),]=sum(postMatrix[c(1,1+length(rounds)),])/(2*length(gridx))}
  postMatrix=sweep(postMatrix,2,colSums(postMatrix),`/`, check.margin = FALSE) #\pi(R_i|X_i,\tau,a,\beta)  
  return(postMatrix)
}

Rprx2=function(postMatrix, rounds, gridx){
  possW=seq(from=plyr::round_any(min(gridx),min(rounds)),to=round_any(max(gridx),min(rounds)),by=min(rounds))
  numbers=matrix(0,nrow=length(rounds)*2,ncol=length(possW),dimnames=list(NULL,possW))
  for(i in 1:length(c(rounds,rounds))){
    atemp=factor(plyr::round_any(gridx,accuracy=c(rounds,rounds)[i],f=round),levels=possW)
    numbers[i,]=tapply(postMatrix[i,],INDEX=list(atemp),FUN=function(x) sum(x,na.rm=T))
  }
  numbers[is.na(numbers)]=0
  sweep(numbers,2,colSums(numbers),`/`, check.margin = FALSE)# int \pi(R_i|W_i) dX_i
}

dbc <- function(gridx, x, bw, boundary=c(0,+Inf)){
  kernels <- sapply(x, function(y) {
    dens <- dnorm(gridx,mean=y,sd=bw)
    dens[gridx<boundary[1]|gridx>boundary[2]]=0
    dens
  })
  kernels <- sweep(kernels,2,colSums(kernels),`/`)      
  rowMeans(kernels)/(gridx[2]-gridx[1])
}
          
rprobsRR2=function(RR,beta,gridx,ranefvalues=0){
            pgivenX=matrix(pnorm(rep(beta*log(abs(gridx))+ranefvalues,length(RR))+rep(RR,each=length(gridx))),
                    nrow=length(RR),byrow=T)
            diff(rbind(0,pgivenX,1))
          }
          
          logLik=function(par, new, rguesstemp, unequal, setBias, rounds,ranefvalues, roundvalues){
  RR=sort(par[1:(length(rounds)-1)])
  Bias=0
  beta=0
  if(setBias==TRUE) {Bias=par[length(rounds)]}
  if(unequal==TRUE) {beta=par[length(par)]}
  pgivenX=Rprx(rounds=rounds,RR=RR, Bias=pnorm(Bias), beta=beta, gridx=unique(new), ranefvalues=ranefvalues,
               roundvalues=roundvalues)+1E-96
  probs <- log(pgivenX[
    cumsum(rep(nrow(pgivenX),length.out=length(new)))-cumsum(((sequence(rle(new)$length)>1)*nrow(pgivenX)))-
      nrow(pgivenX)+rguesstemp])
  return(-sum(probs))
}

logLikRandomGamma=function(ranefvalues, RR, Bias, beta, new, rguesstemp, rounds, tau,roundvalues){
  pgivenX=Rprx(rounds=rounds,RR=RR, Bias=pnorm(Bias), beta=beta, gridx=unique(new), ranefvalues=ranefvalues,roundvalues=roundvalues)
  return(sum(log(pgivenX[cumsum(rep(nrow(pgivenX),length.out=length(new)))-
                           cumsum(((sequence(rle(new)$length)>1)*nrow(pgivenX)))-
                           nrow(pgivenX)+rguesstemp]))+sum(log(dnorm(ranefvalues,mean=0,sd=sqrt(tau)))))
}

#' Kernel density estimation for heaped data
#' @param xheaped heaped values from which to estimate density of x
#' @param rounds rounding values, numeric vector of length >=1
#' @param burnin burn-in sample size
#' @param samples sampling iteration size
#' @param setBias if TRUE a rounding Bias parameter is estimated. For values above 0.5, the respondents
#' are more prone to round down, while for values < 0.5 they are more likely to round up  
#' @param bw bandwidth selector method, defaults to "nrd0" see \code{density} for more options  
#' @param boundary TRUE for positive only data (no positive density for negative values)
#' @param unequal if TRUE a probit model is fitted for the rounding probabilities with log(true value) as regressor
#' @param random if TRUE a random effect probit model is fitted for rounding probabilities
#' @param adjust as in \code{density}, the user can multiply the bandwidth by a certain factor such that bw=adjust*bw
#' @param weights optional numeric vector of sampling weights
#' @param recall if TRUE a recall error is introduced to the heaping model 
#' @param recallParams recall error model parameters expression(nu) and expression(eta). Default is c(1/3, 1/3)
#' @return
#' The function returns a list object with the following objects (besides all input objects):
#' \item{\code{meanPostDensity}}{Vector of Mean Posterior Density}
#' \item{\code{gridx}}{Vector Grid on which density is evaluated}
#' \item{\code{resultDensity}}{Matrix with Estimated Density for each iteration}
#' \item{\code{resultRR}}{Matrix with rounding probability threshold values for each iteration (on probit scale)}
#' \item{\code{resultBias}}{Vector with estimated Bias parameter for each iteration}
#' \item{\code{resultBeta}}{Vector with estimated Beta parameter for each iteration}
#' \item{\code{resultX}}{Matrix of true latent values X estimates}
#' @examples
#' #Simple Rounding   ----------------------------------------------------------
#' xtrue=rnorm(3000)
#' xrounded=round(xtrue)
#' est <- dheaping(xrounded,rounds=1,burnin=20,samples=50)
#' plot(est,trueX=xtrue)
#' 
#' #####################
#' #####Heaping
#' #####################
#' #Real Data Example  ----------------------------------------------------------
#' # Student learning hours per week
#' data(students)
#' xheaped <- as.numeric(na.omit(students$StudyHrs))
#' \dontrun{est <- dheaping(xheaped,rounds=c(1,2,5,10), boundary=TRUE, unequal=TRUE,burnin=20,samples=50)
#' plot(est)
#' summary(est)}
#' 
#' #Simulate Data   ----------------------------------------------------------
#' Sim1 <- createSim.Kernelheaping(n=500, distribution="norm",rounds=c(1,10,100),
#' thresholds=c(-0.5244005, 0.5244005), sd=100)
#' \dontrun{est <- dheaping(Sim1$xheaped,rounds=Sim1$rounds)
#' plot(est,trueX=Sim1$x)}
#' 
#' #Biased rounding
#' Sim2 <- createSim.Kernelheaping(n=500, distribution="gamma",rounds=c(1,2,5,10),
#'                      thresholds=c(-1.2815516, -0.6744898, 0.3853205),downbias=0.2,
#'                      shape=4,scale=8,offset=45)
#' \dontrun{est <- dheaping(Sim2$xheaped, rounds=Sim2$rounds, setBias=T, bw="SJ")
#' plot(est, trueX=Sim2$x)
#' summary(est)
#' tracePlots(est)}
#' 
#' Sim3 <- createSim.Kernelheaping(n=500, distribution="gamma",rounds=c(1,2,5,10),
#' thresholds=c(1.84, 2.64, 3.05), downbias=0.75, Beta=-0.5, shape=4, scale=8)
#' \dontrun{est <- dheaping(Sim3$xheaped,rounds=Sim3$rounds,boundary=TRUE,unequal=TRUE,setBias=T)
#' plot(est,trueX=Sim3$x)}
#' @export 
dheaping <- function(xheaped, rounds, burnin=5, samples=10, setBias=FALSE, weights=NULL,
                     bw= "nrd0", boundary=FALSE, unequal=FALSE, random=FALSE, adjust=1,
                     recall=F, recallParams=c(1/3,1/3)){
  if(is.null(weights)){weights=rep(1/length(xheaped),length(xheaped))}
  weights=weights/sum(weights)
  roundvalues=rounds
  if(min(rounds)==0){rounds[rounds==0]=0.5*min(rounds[rounds>0])}
  weights=weights[order(xheaped)]
  xheaped=xheaped[order(xheaped)] #round to lowest rounding value
  weights=weights[!is.na(xheaped)]
  xheaped=xheaped[!is.na(xheaped)]
  xheapedOriginal <- xheaped
  xheaped=plyr::round_any(xheaped,accuracy=min(rounds)) #round to lowest rounding value
  #Create grid
  stepsize=1/2^(which(diff(range(xheaped))/min(rounds)*2^(1:10)>100)[1]) #at least 100 grid points
  gridx=seq(min(xheaped)-max(rounds)*0.5+stepsize/2*min(rounds),
            max(xheaped)+max(rounds)*0.5-stepsize/2*min(rounds),
            by=min(rounds)*stepsize)
  if(boundary==TRUE|unequal==TRUE){gridx=gridx[gridx>0]}
  #Pilot Estimation
  if(boundary==FALSE){Mestimates <- density(xheaped,from=min(gridx),to=max(gridx),n=length(gridx),bw=max(rounds)*2,
                                            weights=weights)$y}
  if(boundary==TRUE){Mestimates <- dbc(gridx=gridx,x=xheaped,bw=max(rounds)*2)}
  #Starting values for x(new), Rounding values(rguesstemp), und Rounding thresholds(RR), Bias (Bias) and beta
  #rguesstemp<-sapply(1:length(xheaped),function(x) rounds[max(which(xheaped[x]%%rounds==0))])
  rguesstemp=rep(0,length(xheaped))
  heapedvalues=unique(xheaped)
  xtable<-table(xheaped)
  #Starting Values
  new<-xheaped
  RR<-qnorm(1:(length(rounds)-1)/(length(rounds)))
  Bias<-0
  beta<-0
  tau<-0.5
  ranefnames=unique(plyr::round_any(gridx,accuracy=min(rounds)))
  ranef=rnorm(length(ranefnames),sd=sqrt(tau))
  ranefvalues=rep(0,length(gridx))
  
  #Result matrices
  resultDensity=matrix(nrow=samples+burnin,ncol=length(gridx))
  resultX=matrix(nrow=samples+burnin,ncol=length(xheaped))
  resultRR=matrix(nrow=samples+burnin,ncol=length(rounds)-1)
  resultBias=matrix(nrow=samples+burnin,ncol=1)
  resultBeta=matrix(nrow=samples+burnin,ncol=1)
  resultRandom=matrix(nrow=samples+burnin,ncol=length(ranefnames),dimnames=list(1:(samples+burnin),ranefnames))
  resultTau=rep(0,times=samples+burnin)
  
  selectionGrid<-lapply(heapedvalues,function(k){
    selectionupper=lapply(1:length(rounds), function(x) 
      which(gridx>k-rounds[x]*0.5&gridx<=k))
    selectionlower=lapply(1:length(rounds), function(x)
      which(gridx<k+rounds[x]*0.5&gridx>=k))
    selection=c(selectionlower,selectionupper)
  })
  
  for(j in 1:(burnin+samples)){
    Rprs1=Rprx(rounds=rounds,RR=RR, Bias=pnorm(Bias), beta=beta,
               gridx=gridx, ranefvalues=ranefvalues,roundvalues=roundvalues)
    Rprs2=Rprx2(Rprs1, rounds=rounds, gridx=gridx)[,as.character(heapedvalues)]
    if(recall==FALSE){
      getNew <- getNewX(heapedvalues, selectionGrid, Mestimates, Rprs1, Rprs2,
                        xtable, xheaped, gridx, xheapedOriginal, roundvalues, rguesstemp, new)
    }
    if(recall==TRUE){
      getNew <- getNewV(heapedvalues, selectionGrid, Mestimates, Rprs1, Rprs2, xtable, xheaped,
                        gridx, xheapedOriginal, roundvalues, rguesstemp, new, Bias, beta, ranefvalues, RR, recallParams)
    }
    new=getNew[[1]]
    rguesstemp=getNew[[2]]
    #ordered
    newx=getNew[[1]][order(getNew[[1]])]
    rguesstempx=getNew[[2]][order(getNew[[1]])]
    #Laplace-Approximation:
    par=RR
    if(setBias==TRUE){par=c(par,Bias)}
    if(unequal==TRUE){par=c(par,beta)}
    if(length(rounds)>1){
      ranefvaluestemp=ranefvalues[match(unique(newx),gridx)]
      laplace=optim(par,logLik,new=newx,rguesstemp=rguesstempx,
                    unequal=unequal, setBias=setBias,
                    ranefvalues=ranefvaluestemp,rounds=rounds ,hessian=T, roundvalues=roundvalues,
                    method="BFGS")
      par=MASS::mvrnorm(1,mu=laplace$par,Sigma=solve(laplace$hessian+diag(1E-7,length(par))))
      
      #     Metropolis-Hastings als Alternative
      #       parNew=par+rnorm(length(par),sd=0.05)
      #       alpha=min(1,exp(-logLik(par = parNew, new=new,rguesstemp=rguesstemp, unequal=unequal, setBias=setBias,
      #                              postMatrix=matrix(0,ncol=length(unique(new)),nrow=2*length(rounds),dimnames=list(NULL,unique(new[order(new)]))),
      #                              ranefvalues=ranefvaluestemp,rounds=rounds)+
      #                         logLik(par = par, new=new,rguesstemp=rguesstemp, unequal=unequal, setBias=setBias,
      #                                postMatrix=matrix(0,ncol=length(unique(new)),nrow=2*length(rounds),dimnames=list(NULL,unique(new[order(new)]))),
      #                                ranefvalues=ranefvaluestemp,rounds=rounds)
      #       ))
      #       if(runif(1)<alpha){par=parNew}
      
      if(setBias==TRUE) {Bias=par[length(rounds)]}
      if(unequal==TRUE) {beta=par[length(par)]}
      RR=sort(par[1:(length(rounds)-1)])
      if(random==TRUE){
        #Metropolis Hastings for RE
        ranefnamestemp=unique(plyr::round_any(newx,accuracy=min(rounds)))
        raneftemp=ranef[match(ranefnamestemp,ranefnames)]
        ranefnew=raneftemp
        for(k in 1:length(raneftemp)){
          ranefnew[k]=raneftemp[k]+rnorm(1,sd=0.2)
          alpha=min(1,exp(logLikRandomGamma(ranefvalues=ranefnew[k], RR=RR, Bias=Bias, beta=beta,
                                            new=newx[plyr::round_any(newx,accuracy=min(rounds))==ranefnamestemp[k]],
                                            rguesstemp=rguesstempx[plyr::round_any(newx,accuracy=min(rounds))==ranefnamestemp[k]],
                                            rounds=rounds, tau=tau, roundvalues=roundvalues)-
                            logLikRandomGamma(ranefvalues=raneftemp[k], RR=RR, Bias=Bias, beta=beta,
                                              new=newx[plyr::round_any(newx,accuracy=min(rounds))==ranefnamestemp[k]],
                                              rguesstemp=rguesstempx[plyr::round_any(newx,accuracy=min(rounds))==ranefnamestemp[k]],
                                              rounds=rounds, tau=tau ,roundvalues=roundvalues)
          ))
          if(is.na(alpha)){alpha=0}
          if(runif(1)<alpha){raneftemp[k]=ranefnew[k]}
        }
        tau=1/rgamma(1,shape=0.001+length(raneftemp)/2,scale=(0.001+0.5*sum((raneftemp)^2))^-1)
        ranef[match(ranefnamestemp,ranefnames)]=raneftemp
        ranefvalues=ranef[match(plyr::round_any(gridx,accuracy=min(rounds)),ranefnames)]
      }
    }
    if(0 %in% roundvalues){
      new[xheapedOriginal%%min(roundvalues[roundvalues>0])!=0]=xheapedOriginal[xheapedOriginal%%min(roundvalues[roundvalues>0])!=0]
    }
    if(recall==T){
      sdV <- recallParams[1]*c(rounds,rounds)[rguesstemp]+
        recallParams[2]*sapply(new,function(x) 
          diff(c(0,pnorm(rep(beta*log(abs(x)),length(RR))+RR),1))%*%
            roundvalues)
      
      new=sapply(1:length(new), function(y) sample(x=gridx,size=1,
                                                   prob=dnorm(gridx,mean=new[y]+qnorm(pnorm(Bias))*sdV[y],
                                                              sd=sdV[y])*Mestimates))
    }

    h <- density(new,bw=bw)$bw*adjust
    if(boundary==TRUE){Mestimates <- dbc(gridx=gridx,x=new,bw=h)}
    if(boundary==FALSE){Mestimates <- density(new,from=min(gridx),to=max(gridx),n=length(gridx),bw=h,weights=weights)$y}
    #Save Results
    resultDensity[j,]=Mestimates
    resultRR[j,]=RR
    resultBias[j,]=pnorm(Bias)
    resultBeta[j,]=beta
    resultX[j,]=new
    resultTau[j]=tau
    resultRandom[j,]=ranef
    print(paste("Iteration:",j,"of", burnin+samples))
  }
  meanPostDensity=apply(resultDensity[-c(1:burnin),],2,mean)
  est<-list(meanPostDensity=meanPostDensity,resultDensity=resultDensity,resultRR=resultRR, resultBias=resultBias,
            resultBeta=resultBeta, resultX=resultX, xheaped=xheaped, gridx=gridx, boundary=boundary, rounds=roundvalues,
            setBias=setBias, unequal=unequal, bw=bw, burnin=burnin, samples=samples, adjust=adjust, 
            resultTau=resultTau, resultRandom=resultRandom, random=random, weights=weights)
  class(est) <- "Kernelheaping"
  return(est)
}

#' Plot Kernel density estimate of heaped data naively and corrected by partly bayesian model
#' @param x Kernelheaping object produced by \code{dheaping} function
#' @param trueX optional, if true values X are known (in simulations, for example) the 'Oracle' density estimate is added as well
#' @param ... additional arguments given to standard plot function
#' @return plot with Kernel density estimates (Naive, Corrected and True (if provided))
#' @method plot Kernelheaping
#' @export
plot.Kernelheaping <- function(x,trueX=NULL, ...){
  if(x$boundary==FALSE){
    plot(density(x$xheaped,bw=x$bw,adjust=x$adjust, weights=x$weights),xlab="x",ylab="Density",main="", ...)
    if(!is.null(trueX)){lines(density(trueX,bw=x$bw,adjust=x$adjust, weights=x$weights),col="blue",lty=2,lwd=2)}
  }
  if(x$boundary==TRUE){
    plot(dbc(gridx=x$gridx,x=x$xheaped,bw=density(x$xheaped,bw=x$bw,adjust=x$adjust)$bw)~x$gridx,type="l",xlab="x",ylab="Density", main="", ...)
    if(!is.null(trueX)){lines(dbc(x$gridx, trueX, bw=density(trueX,bw=x$bw,adjust=x$adjust)$bw)~x$gridx,col="blue",lty=2,lwd=2)}
  }
  lines(x$meanPostDensity~x$gridx,col="red",lty=4,lwd=2)
  if(is.null(trueX)){legend("topright",c("Naive","Corrected"),col=c("black","red"),lty=c(1,4,2))}
  if(!is.null(trueX)){legend("topright",c("Naive","Corrected","Oracle"),col=c("black","red","blue"),lty=c(1,4,2),lwd=c(1,2,2))}
}

#' Prints some descriptive statistics (means and quantiles) for the estimated rounding, bias and acceleration (beta) parameters
#' @param object Kernelheaping object produced by \code{dheaping} function
#' @param ... unused
#' @return Prints summary statistics
#' @method summary Kernelheaping
#' @export
summary.Kernelheaping <- function(object, ...){
  Threshold=c(-Inf,colMeans(object$resultRR[-c(1:object$burnin),]))
  names(Threshold)=object$rounds
  cat("Mean Posterior Threshold Values:", "\n")
  cat(Threshold, "\n")
  if(object$setBias==TRUE){
    Bias=c(mean(object$resultBias[-c(1:object$burnin)]),quantile(object$resultBias[-c(1:object$burnin)],c(0.025,0.975)))
    names(Bias)[1]="Mean"
    cat("\n", "Bias, Mean and 95% credible interval:", "\n")
    cat(Bias, "\n")
  }
  if(object$unequal==TRUE){
    unequal=c(mean(object$resultBeta[-c(1:object$burnin)]),quantile(object$resultBeta[-c(1:object$burnin)],c(0.025,0.975)))
    names(unequal)[1]="Mean"
    cat("\n", "Beta, Mean and 95% credible interval:", "\n")
    cat(unequal, "\n")}
}

#' Plots some trace plots for the rounding, bias and acceleration (beta) parameters
#' @param x Kernelheaping object produced by \code{dheaping} function
#' @param ... additional arguments given to standard plot function
#' @return Prints summary statistics
#' @export
tracePlots <- function(x, ...){
  par(mfrow=c(2,ceiling(ncol(x$resultRR)/2)))
  for(i in 1:ncol(x$resultRR)){
    plot(x$resultRR[,i],xlab="iteration", ylab="value", main=paste("Traceplot Threshold", i), type="l")
  }
  par(mfrow=c(1,1))
  if(x$setBias==TRUE){
    plot(x$resultBias,xlab="iteration", ylab="value", main="Traceplot Bias", type="l")
  }
  
  if(x$unequal==TRUE){
    plot(x$resultBeta,xlab="iteration", ylab="value", main="Traceplot Beta", type="l")
  }
  if(x$random==TRUE){
    plot(x$resultTau,xlab="iteration", ylab="value", main="Traceplot tau", type="l")
  }
}

#' Simulation of heaping correction method
#' @param simRuns number of simulations runs
#' @param n sample size
#' @param distribution name of the distribution where random sampling is available, e.g. "norm"
#' @param rounds rounding values, numeric vector of length >=1
#' @param thresholds rounding thresholds
#' @param burnin burn-in sample size
#' @param samples sampling iteration size
#' @param downbias Bias parameter used in the simulation
#' @param setBias if TRUE a rounding Bias parameter is estimated. For values above 0.5, the respondents
#' are more prone to round down, while for values < 0.5 they are more likely to round up  
#' @param bw bandwidth selector method, defaults to "nrd0" see \code{density} for more options  
#' @param offset location shift parameter used simulation in simulation
#' @param boundary TRUE for positive only data (no positive density for negative values)
#' @param Beta Parameter of the probit model for rounding probabilities used in simulation
#' @param unequal if TRUE a probit model is fitted for the rounding probabilities with log(true value) as regressor
#' @param adjust as in \code{density}, the user can multiply the bandwidth by a certain factor such that bw=adjust*bw
#' @param ... additional attributes handed over to \code{createSim.Kernelheaping}
#' @return List of estimation results
#' @examples
#' \dontrun{Sims1 <- sim.Kernelheaping(simRuns=2, n=500, distribution="norm", 
#' rounds=c(1,10,100), thresholds=c(0.3,0.4,0.3), sd=100)}
#' @export
sim.Kernelheaping <- function(simRuns, n, distribution, rounds, thresholds, downbias=0.5, setBias = FALSE, Beta=0, unequal = FALSE, burnin = 5, samples = 10,
                              bw = "nrd0", offset=0, boundary = FALSE, adjust = 1, ...){
  lapply(1:simRuns, function(x){
    print(x)
    Sim <- createSim.Kernelheaping(n, distribution, rounds, thresholds, downbias=downbias, Beta=Beta, offset=offset, ...)
    est <- dheaping(Sim$xheaped, rounds=Sim$rounds, burnin = burnin, samples = samples, setBias = setBias,
                    bw = bw, boundary = boundary, unequal = unequal, adjust = adjust)
    
    if(est$boundary==T){
      naiveD  <- dbc(est$gridx, est$xheaped, density(est$xheaped,from=min(est$gridx),to=max(est$gridx),n=length(est$gridx),bw=est$bw,adjust=est$adjust)$bw)
      oracleD <- dbc(est$gridx, Sim$x, density(Sim$x,from=min(est$gridx),to=max(est$gridx),n=length(est$gridx),bw=est$bw,adjust=est$adjust)$bw)
    }
    if(est$boundary==F){
      naiveD  <- density(est$xheaped,from=min(est$gridx),to=max(est$gridx),n=length(est$gridx),bw=est$bw,adjust=est$adjust)$y
      oracleD <- density(Sim$x,from=min(est$gridx),to=max(est$gridx),n=length(est$gridx),bw=est$bw,adjust=est$adjust)$y
    }
    trueD <- do.call(paste("d",distribution,sep=""), list(est$grid-offset, ...))
    return(list(est=est,Sim=Sim,naiveD=naiveD,oracleD=oracleD,trueD=trueD))
  })
}
#' Simulation Summary
#' @param sim Simulation object returned from sim.Kernelheaping
#' @param coverage probability for computing coverage intervals
#' @return list with summary statistics
#' @export
simSummary.Kernelheaping <- function(sim, coverage=0.9){
  RMISE <- sqrt(sapply(sim,function(x) cbind(sum(diff(x$est$grid)[1]*(x$est$meanPostDensity-x$trueD)^2),
                                             sum(diff(x$est$grid)[1]*(x$naiveD-x$trueD)^2),
                                             sum(diff(x$est$grid)[1]*(x$oracleD-x$trueD)^2))))
  RRmeans <- sapply(sim,function(x) colMeans(as.matrix(x$est$resultRR[-c(1:x$est$burnin),],nrow=x$est$samples)))
  RRlower <- sapply(sim,function(x) apply(x$est$resultRR[-c(1:x$est$burnin),],2,quantile,0.5*(1-coverage)))
  RRupper <- sapply(sim,function(x) apply(x$est$resultRR[-c(1:x$est$burnin),],2,quantile,1-(0.5*(1-coverage))))
  Biasmean <- sapply(sim,function(x) mean(x$est$resultBias[-c(1:x$est$burnin),]))
  Biaslower <- sapply(sim,function(x) quantile(x$est$resultBias[-c(1:x$est$burnin),],0.5*(1-coverage)))
  Biasupper <- sapply(sim,function(x) quantile(x$est$resultBias[-c(1:x$est$burnin),],1-(0.5*(1-coverage))))
  Betamean <- sapply(sim,function(x) mean(x$est$resultBeta[-c(1:x$est$burnin),]))
  Betalower <- sapply(sim,function(x) quantile(x$est$resultBeta[-c(1:x$est$burnin),],0.5*(1-coverage)))
  Betaupper <- sapply(sim,function(x) quantile(x$est$resultBeta[-c(1:x$est$burnin),],1-(0.5*(1-coverage))))
  #coverage 
  coverageRR <- rowSums(sapply(1:length(sim), function(x) sim[[x]]$Sim$thresholds > RRlower[,x] & sim[[x]]$Sim$thresholds < RRupper[,x]))/length(sim)
  coverageBias <- sum(sapply(1:length(sim), function(x) sim[[x]]$Sim$downbias>Biaslower[x] & sim[[x]]$Sim$downbias<Biasupper[x]))/length(sim)
  coverageBeta <- sum(sapply(1:length(sim), function(x) sim[[x]]$Sim$Beta>Betalower[x] & sim[[x]]$Sim$Beta<Betaupper[x]))/length(sim)
  return(list(RMISE=RMISE,RRmeans=RRmeans,RRlower=RRlower,RRupper=RRupper,coverageRR=coverageRR,
              Biasmean=Biasmean,Biaslower=Biaslower,Biasupper=Biasupper,coverageBias=coverageBias,
              Betamean=Betamean,Betalower=Betalower,Betaupper=Betaupper,coverageBeta=coverageBeta))
}
#' Create heaped data for Simulation
#' @param n sample size
#' @param distribution name of the distribution where random sampling is available, e.g. "norm"
#' @param rounds rounding values
#' @param thresholds rounding thresholds (for Beta=0)
#' @param offset certain value added to all observed random samples
#' @param downbias bias parameter
#' @param Beta acceleration paramter
#' @param ... additional attributes handed over to "rdistribution" (i.e. rnorm, rgamma,..)
#' @return List of heaped values, true values and input parameters
#' @export
createSim.Kernelheaping <- function(n, distribution, rounds, thresholds, offset=0, downbias=0.5, Beta=0, ...){
  xtrue=do.call(paste("r",distribution,sep=""), list(n, ...))+offset
  RR0=thresholds
  down=lapply(xtrue,function(x) (x%%rounds<=rounds*0.5)*downbias*rprobsRR2(RR=RR0,beta=Beta,gridx=x)+
                (x%%rounds>=rounds*0.5)*(1-downbias)*rprobsRR2(RR=RR0,beta=Beta,gridx=x))
  down=lapply(down, function(x) x/sum(x))
  roundings=sapply(down,function(z) sample(rounds,size=1,replace=TRUE,prob=z))
  xheaped=plyr::round_any(xtrue,accuracy=roundings, f=round)
  return(list(xheaped=xheaped, rounds=rounds, thresholds=thresholds, downbias=downbias, Beta=Beta, x=xtrue))
}
#' Bivariate kernel density estimation for rounded data
#' @param xrounded rounded values from which to estimate bivariate density, matrix with 2 columns (x,y) 
#' @param roundvalue rounding value (side length of square in that the true value lies around the rounded one)
#' @param burnin burn-in sample size
#' @param samples sampling iteration size
#' @param adaptive set to TRUE for adaptive bandwidth
#' @param gridsize number of evaluation grid points
#' @return
#' The function returns a list object with the following objects (besides all input objects):
#' \item{\code{Mestimates}}{kde object containing the corrected density estimate}
#' \item{\code{gridx}}{Vector Grid on which density is evaluated (x)}
#' \item{\code{gridy}}{Vector Grid on which density is evaluated (y)}
#' \item{\code{resultDensity}}{Array with Estimated Density for each iteration}
#' \item{\code{resultX}}{Matrix of true latent values X estimates}
#' \item{\code{delaigle}}{Matrix of Delaigle estimator estimates}
#' @examples
#' # Create Mu and Sigma  -----------------------------------------------------------
#' mu1 <- c(0, 0)
#' mu2 <- c(5, 3)
#' mu3 <- c(-4, 1)
#' Sigma1 <- matrix(c(4, 3, 3, 4), 2, 2)
#' Sigma2 <- matrix(c(3, 0.5, 0.5, 1), 2, 2)
#' Sigma3 <- matrix(c(5, 4, 4, 6), 2, 2)
#' # Mixed Normal Distribution -------------------------------------------------------
#' mus <- rbind(mu1, mu2, mu3)
#' Sigmas <- rbind(Sigma1, Sigma2, Sigma3)
#' props <- c(1/3, 1/3, 1/3)
#' \dontrun{xtrue=rmvnorm.mixt(n=1000, mus=mus, Sigmas=Sigmas, props=props)
#' roundvalue=2
#' xrounded=plyr::round_any(xtrue,roundvalue)
#' est <- dbivr(xrounded,roundvalue=roundvalue,burnin=5,samples=10)
#' 
#' #Plot corrected and Naive distribution
#' plot(est,trueX=xtrue)
#' #for comparison: plot true density
#'  dens=dmvnorm.mixt(x=expand.grid(est$Mestimates$eval.points[[1]],est$Mestimates$eval.points[[2]]),
#'   mus=mus, Sigmas=Sigmas, props=props)
#'  dens=matrix(dens,nrow=length(est$gridx),ncol=length(est$gridy))
#'  contour(dens,x=est$Mestimates$eval.points[[1]],y=est$Mestimates$eval.points[[2]],
#'     xlim=c(min(est$gridx),max(est$gridx)),ylim=c(min(est$gridy),max(est$gridy)),main="True Density")}
#'@export
dbivr <- function(xrounded, roundvalue, burnin=2, samples=5, adaptive=FALSE, gridsize=1000){
  #Create grid - sparr package requires identical grid length on both axis:
  gridx=seq(min(xrounded[,1])-0.5*roundvalue,max(xrounded[,1])+0.5*roundvalue,length=150)
  gridy=seq(min(xrounded[,2])-0.5*roundvalue,max(xrounded[,2])+0.5*roundvalue,length=150)
  
  #Pilot Estimation
  Mestimates <- ks::kde(x=xrounded, H=diag(c(roundvalue,roundvalue))^2,gridsize=c(length(gridx),length(gridy)),
                        xmin=c(min(gridx),min(gridy)),xmax=c(max(gridx),max(gridy)))
  
  #Result matrices
  resultDensity=array(dim=c(burnin+samples,length(gridx),length(gridy)))
  resultX=array(dim=c(samples+burnin,nrow(xrounded),2))
  
  rvalues=unique(xrounded) #unique rounded values in data
  selectionGrid<-lapply(1:nrow(rvalues),function(k){
    selectionX=which(Mestimates$eval.points[[1]]>=rvalues[k,1]-roundvalue*0.5&Mestimates$eval.points[[1]]<rvalues[k,1]+roundvalue*0.5)
    selectionY=which(Mestimates$eval.points[[2]]>=rvalues[k,2]-roundvalue*0.5&Mestimates$eval.points[[2]]<rvalues[k,2]+roundvalue*0.5)
    list(selectionX,selectionY)
  })
  
  
  #Delaigle Estimator
  delaigle=aggregate(list(length=rep(1,nrow(xrounded))), data.frame(xrounded), length)
  delaigle[,3]=delaigle[,3]/sum(delaigle[,3])/roundvalue^2
  delaigleest=matrix(0,nrow=length(gridx),ncol=length(gridy))
  for(i in 1:nrow(delaigle)){
    x=which(delaigle[i,1]==rvalues[,1]&delaigle[i,2]==rvalues[,2])
    delaigleest[selectionGrid[[x]][[1]],selectionGrid[[x]][[2]]]=delaigle[i,3]
  }
  
  #SEM Estimator
  for(j in 1:(burnin+samples)){
    new=c()
    for(i in 1:nrow(rvalues)){
      probs=as.vector(Mestimates$estimate[selectionGrid[[i]][[1]],selectionGrid[[i]][[2]]])
      points=cbind(rep(Mestimates$eval.points[[1]][selectionGrid[[i]][[1]]],
                       times=length(Mestimates$eval.points[[2]][selectionGrid[[i]][[2]]])),
                   rep(Mestimates$eval.points[[2]][selectionGrid[[i]][[2]]],
                       each=length(Mestimates$eval.points[[1]][selectionGrid[[i]][[1]]])))    
      npoints=length(which(xrounded[,1]==rvalues[i,1]&xrounded[,2]==rvalues[i,2]))
      new=rbind(new,points[sample(1:nrow(points),size=npoints,replace=T,prob=probs),])
    }
    
    #recompute H
    if(adaptive==FALSE){
      H <- ks::Hpi(x = new, binned=TRUE)
    }
    if(adaptive==TRUE){
      H <- ks::Hpi(x = new, binned=TRUE)
      H <- sqrt(sqrt(H[1,1]*H[2,2]))
    }
    #recompute density
    if(adaptive==FALSE){
      Mestimates <- ks::kde(x=new, H=H,gridsize=c(length(gridx),length(gridy)),bgridsize=c(length(gridx),length(gridy)),
                            xmin=c(min(gridx),min(gridy)),xmax=c(max(gridx),max(gridy)),binned=TRUE)
    }
    if(adaptive==TRUE){
      counts <- plyr::count(new)
      MestimatesAd <- sparr::bivariate.density(data=counts[,c(1:2)],pilotH=H,res=length(gridx),xrange=range(gridx),
                                               yrange=range(gridy),adaptive=TRUE,comment=FALSE,counts=counts[,3])
      Mestimates$estimate=MestimatesAd$Zm
    }
    Mestimates$estimate[is.na(Mestimates$estimate)]=1E-96
    resultDensity[j,,]=Mestimates$estimate
    resultX[j,,]=new
    print(paste("Iteration:",j,"of", burnin+samples))
  }
  Mestimates$estimate=apply(resultDensity[-c(1:burnin),,],c(2,3),mean)
  est<-list(Mestimates=Mestimates,resultDensity=resultDensity,resultX=resultX,
            xrounded=xrounded, gridx=gridx, gridy=gridy, roundvalue=roundvalue,
            burnin=burnin, samples=samples, adaptive=adaptive, delaigle=delaigleest)
  class(est) <- "bivrounding"
  return(est)
}

#' Plot Kernel density estimate of heaped data naively and corrected by partly bayesian model
#' @param x bivrounding object produced by \code{dbivr} function
#' @param trueX optional, if true values X are known (in simulations, for example) the 'Oracle' density estimate is added as well
#' @param ... additional arguments given to standard plot function
#' @return plot with Kernel density estimates (Naive, Corrected and True (if provided))
#' @method plot bivrounding
#' @export
plot.bivrounding <- function(x, trueX=NULL, ...){
  par(mfrow=c(2,2))
  
  #select Levels
  if(x$adaptive==FALSE){
    Naive=kde(x$xrounded,gridsize=c(length(x$gridx),length(x$gridy)),
              xmin=c(min(x$gridx),min(x$gridy)),xmax=c(max(x$gridx),max(x$gridy)))
    contour(x$Mestimates$estimate,x=x$Mestimates$eval.points[[1]],y=x$Mestimates$eval.points[[2]],
            col=c("white",rainbow(6,start=0.5,end=0.6)),main="Corrected")
    contour(Naive$estimate,x=x$Mestimates$eval.points[[1]],y=x$Mestimates$eval.points[[2]],
            col=c("white",rainbow(6,start=0.5,end=0.6)),main="Naive")
    if(!is.null(trueX)){
      Oracle=kde(trueX,gridsize=c(length(x$gridx),length(x$gridy)),
                 xmin=c(min(x$gridx),min(x$gridy)),xmax=c(max(x$gridx),max(x$gridy)))
      contour(Oracle$estimate,x=x$Mestimates$eval.points[[1]],y=x$Mestimates$eval.points[[2]],
              col=c("white",rainbow(6,start=0.5,end=0.6)),main="Oracle")
    }
  }
  
  if(x$adaptive==TRUE){
    Naive=sparr::bivariate.density(data=x$xrounded,pilotH=ks::hpi(x = x$xrounded,binned=TRUE),res=length(x$gridx),xrange=range(x$gridx),
                                   yrange=range(x$gridy),adaptive=TRUE,comment=FALSE)
    contour(x$Mestimates$estimate,x=x$Mestimates$eval.points[[1]],y=x$Mestimates$eval.points[[2]],
            col=c("white",rainbow(6,start=0.5,end=0.6)),main="Corrected")
    contour(Naive$Zm,x=x$Mestimates$eval.points[[1]],y=x$Mestimates$eval.points[[2]],
            col=c("white",rainbow(6,start=0.5,end=0.6)),main="Naive")
    if(!is.null(trueX)){
      Oracle=sparr::bivariate.density(data=trueX,pilotH=ks::hpi(x = trueX),res=length(x$gridx),xrange=range(x$gridx),
                                      yrange=range(x$gridy),adaptive=TRUE,comment=FALSE)
      contour(Oracle$Zm,x=x$Mestimates$eval.points[[1]],y=x$Mestimates$eval.points[[2]],
              col=c("white",rainbow(6,start=0.5,end=0.6)),main="Oracle")
    }
  }
  
  image(x$delaigle,oldstyle=TRUE,x=x$Mestimates$eval.points[[1]],xlab="",ylab="",
        y=x$Mestimates$eval.points[[2]],col=c("white",rainbow(6,start=0.5,end=0.6)),main="Delaigle")
  box()
  par(mfrow=c(1,1))
}

getNewX <- function(heapedvalues, selectionGrid, Mestimates, Rprs1, Rprs2,
                    xtable, xheaped, gridx, xheapedOriginal, roundvalues, rguesstemp, new){
  for(i in 1:length(heapedvalues)){
    selection=selectionGrid[[i]]
    selectionprobs=lapply(1:length(selection),function(x)
      Mestimates[selection[[x]]]*
        #f(X_i,R_i|.)
        Rprs1[x,selection[[x]]]/(sum(Rprs1[x,selection[[x]]]))*
        (Rprs2[x,i]))
    selectionprobs <- unlist(selectionprobs)
    selectionprobs[is.na(selectionprobs)] <- 1e-16
    temprounds=unlist(lapply(1:length(selection), function(x) rep(x,times=length(selection[[x]]))))
    temp=sample(1:length(selectionprobs),size=xtable[i],prob=selectionprobs,replace=TRUE)
    rguesstemp[xheaped==heapedvalues[i]]=temprounds[temp]
    new[xheaped==heapedvalues[i]]=gridx[unlist(selection)[temp]]
  }
  if(0 %in% roundvalues){
    rguesstemp[xheapedOriginal%%min(roundvalues[roundvalues>0])!=0]=sample(c(1,length(roundvalues)+1),
                                                                           length(rguesstemp[xheapedOriginal%%min(roundvalues[roundvalues>0])!=0]),replace=T)
    rguesstemp[xheapedOriginal%%min(roundvalues[roundvalues>0])==0&rguesstemp %in% c(1,length(roundvalues)+1)] <- 
      rguesstemp[xheapedOriginal%%min(roundvalues[roundvalues>0])==0&rguesstemp %in% c(1,length(roundvalues)+1)]+1
  }
  return(list(new,rguesstemp))
}

getNewV <- function(heapedvalues, selectionGrid, Mestimates, Rprs1, Rprs2, xtable, xheaped,
                    gridx, xheapedOriginal, roundvalues, rguesstemp, new, Bias, beta, ranefvalues, RR, recallParams,rounds){
  for(i in 1:length(heapedvalues)){
    selection=selectionGrid[[i]]
    selectionprobs=lapply(1:length(selection),function(x)
      #f(X_i,R_i|.)
      Rprs1[x,selection[[x]]]/(sum(Rprs1[x,selection[[x]]]))*
        (Rprs2[x,i]))
    selectionprobs <- unlist(selectionprobs)
    selectionprobs[is.na(selectionprobs)] <- 1e-16
    
    temprounds=unlist(lapply(1:length(selection), function(x) rep(x,times=length(selection[[x]]))))
    
    sdV <- recallParams[1]*c(roundvalues,roundvalues)[temprounds]+
      recallParams[2]*sapply(unlist(selection),function(x) 
        diff(c(0,pnorm(rep(beta*log(abs(gridx[x])),length(RR))+RR),1))%*%roundvalues)
    
    
    selectionprobs2=sapply(new[xheaped==heapedvalues[i]],function(x)
      dnorm(gridx[unlist(selection)],mean=x-qnorm(pnorm(Bias))*sdV,
            sd=sdV))
    
    temp=sapply(1:ncol(selectionprobs2), function(j)
      sample(1:length(selectionprobs),size=1,prob=selectionprobs*selectionprobs2[,j],
             replace=TRUE))
    
    rguesstemp[xheaped==heapedvalues[i]]=temprounds[temp]
    new[xheaped==heapedvalues[i]]=gridx[unlist(selection)[temp]]
  }
  if(0 %in% roundvalues){
    rguesstemp[xheapedOriginal%%min(roundvalues[roundvalues>0])!=0]=sample(c(1,length(roundvalues)+1),
                                                                           length(rguesstemp[xheapedOriginal%%min(roundvalues[roundvalues>0])!=0]),replace=T)
    rguesstemp[xheapedOriginal%%min(roundvalues[roundvalues>0])==0&rguesstemp %in% c(1,length(roundvalues)+1)] <- 
      rguesstemp[xheapedOriginal%%min(roundvalues[roundvalues>0])==0&rguesstemp %in% c(1,length(roundvalues)+1)]+1
  }
  return(list(new,rguesstemp))
}

#' Kernel density estimation for classified data
#' @param xclass classified values; factor with ordered factor values
#' @param classes numeric vector of classes; Inf as last value is allowed
#' @param burnin burn-in sample size
#' @param samples sampling iteration size
#' @param boundary TRUE for positive only data (no positive density for negative values)
#' @param bw bandwidth selector method, defaults to "nrd0" see \code{density} for more options  
#' @param evalpoints number of evaluation grid points
#' @param adjust as in \code{density}, the user can multiply the bandwidth by a certain factor such that bw=adjust*bw
#' @return
#' The function returns a list object with the following objects (besides all input objects):
#' \item{\code{Mestimates}}{kde object containing the corrected density estimate}
#' \item{\code{gridx}}{Vector Grid on which density is evaluated}
#' \item{\code{resultDensity}}{Matrix with Estimated Density for each iteration}
#' \item{\code{resultX}}{Matrix of true latent values X estimates}
#' @examples
#' x=rlnorm(500, meanlog = 8, sdlog = 1)
#' classes <- c(0,500,1000,1500,2000,2500,3000,4000,5000,6000,8000,10000,15000,Inf)
#' xclass <- cut(x,breaks=classes)
#' densityEst <- dclass(xclass=xclass, classes=classes, burnin=2, samples=5, evalpoints=1000)
#' hist(densityEst$xclass,breaks=densityEst$classes)
#' lines(densityEst$Mestimates~densityEst$gridx,col="purple",lwd=2)
#' @export 
dclass <- function(xclass, classes, burnin=2, samples=5, boundary=FALSE, bw="nrd0",
                   evalpoints=200, adjust=1){
  if(max(classes)==Inf){classes[length(classes)]=3*classes[length(classes)-1]}
  classmeans <- sapply(1:(length(classes)-1), function(x) 1/2*(classes[x+1]+classes[x]))
  levels(xclass) <- classmeans
  lengths=as.vector(table(xclass))
  
  xclass <- as.numeric(as.character(xclass))
  #Create grid - sparr package requires identical grid length on both axis:
  gridx=seq(min(classes),max(classes),length=evalpoints)
  #Pilot Estimation
  if(boundary==FALSE){Mestimates <- density(xclass,from=min(gridx),to=max(gridx),n=length(gridx),bw=2*max(classes)/length(classes))$y}
  if(boundary==TRUE){Mestimates <- dbc(gridx=gridx,x=xclass,bw=2*max(classes)/length(classes))}
  
  #Result matrices
  resultDensity=matrix(ncol=c(burnin+samples),nrow=length(gridx))
  resultX=matrix(ncol=c(burnin+samples),nrow=length(xclass))
  
  selectionGrid<-lapply(1:(length(classes)-1),function(k){
    selection=which(gridx>=classes[k]&gridx<classes[k+1])
    selection})
  #SEM Estimator
  for(j in 1:(burnin+samples)){
    new=c()
    for(i in 1:(length(classes)-1)){
      probs=as.vector(Mestimates[selectionGrid[[i]]])
      points=gridx[selectionGrid[[i]]]
      npoints=lengths[i]
      new=c(new,points[sample(1:length(points),size=npoints,replace=T,prob=probs)])
    }
    NewDensity <- density(new,from=min(gridx),to=max(gridx),n=length(gridx),bw=bw,adjust=adjust)
    Mestimates <- NewDensity$y
    if(boundary==TRUE){Mestimates <- dbc(gridx=gridx,x=new,bw=NewDensity$bw)}
    resultDensity[,j]=Mestimates
    resultX[,j]=new
    print(paste("Iteration:",j,"of", burnin+samples))
  }
  Mestimates=apply(resultDensity[,-c(1:burnin)],c(1),mean)
  est<-list(Mestimates=Mestimates,resultDensity=resultDensity,resultX=resultX,
            xclass=xclass, gridx=gridx, classes=classes,
            burnin=burnin, samples=samples)
  class(est) <- "classDensity"
  return(est)
}

#' Bivariate Kernel density estimation for data classified in polygons or shapes
#' @param data matrix with at least 3 columns: x-coordinate, y-coordinate (i.e. center of polygon) and number of observations in area Optional fourth column: ID-Variable if area consists of more than 1 polygon
#' @param burnin burn-in sample size
#' @param samples sampling iteration size
#' @param adaptive TRUE for adaptive kernel density estimation
#' @param shapefile shapefile with number of polygons equal to nrow(data)
#' @param gridsize number of evaluation grid points
#' @return
#' The function returns a list object with the following objects (besides all input objects):
#' \item{\code{Mestimates}}{kde object containing the corrected density estimate}
#' \item{\code{gridx}}{Vector Grid of x-coordinates on which density is evaluated}
#' \item{\code{gridy}}{Vector Grid of y-coordinates on which density is evaluated}
#' \item{\code{resultDensity}}{Matrix with Estimated Density for each iteration}
#' \item{\code{resultX}}{Matrix of true latent values X estimates}
#' @export 
dshapebivr <- function(data, burnin=2, samples=5, adaptive=FALSE, shapefile, gridsize=200){
  ###########################################################
  ##### get polygon shape coordinates from data #####
  ###########################################################
  pol.x <- list()
  pol.y <- list()
  for (i in 1:length(shapefile@polygons)) {
    pol.x[[i]] <- shapefile@polygons[[i]]@Polygons[[1]]@coords[,1]
    pol.y[[i]] <- shapefile@polygons[[i]]@Polygons[[1]]@coords[,2]
  }
  if(ncol(data)==3){           # Create ID Variable if not present
    data[,4] <- 1:nrow(data)
  }
  
  npoints <- data[!(duplicated(data[,4])),3] #delete duplicated points
  
  #Create grid - sparr package requires identical grid length on both axis:
  gridx=seq(shapefile@bbox[1,1],shapefile@bbox[1,2],length=gridsize)
  gridy=seq(shapefile@bbox[2,1],shapefile@bbox[2,2],length=gridsize)
  grid <- as.matrix(expand.grid(gridx,gridy))
  
  #Pilot Estimation
  Mestimates <- ks::kde(x=data[,c(1,2)], H=diag(c(diff(shapefile@bbox[1,])/sqrt(length(shapefile@polygons)),
                                                  c(diff(shapefile@bbox[2,])/sqrt(length(shapefile@polygons)))))^2,
                        gridsize=c(gridsize,gridsize),
                        xmin=c(min(gridx),min(gridy)),
                        xmax=c(max(gridx),max(gridy)),w=nrow(data)*data[,3]/sum(data[,3]))
  #Result matrices
  resultDensity=array(dim=c(burnin+samples,length(gridx),length(gridy)))
  resultX=array(dim=c(samples+burnin,sum(npoints),2))
  
  selectionGrid <- lapply(1:length(unique(data[,4])),function(i){
    y=which(data[,4]==unique(data[,4])[i])
    as.vector(unlist(sapply(y,function(y) which(point.in.polygon(grid[,1],grid[,2],pol.x[[y]],pol.y[[y]])==1))))
  })
  
  #SEM Estimator
  for(j in 1:(burnin+samples)){
    new=c()
    for(i in 1:length(selectionGrid)){
      probs=Mestimates$estimate[cbind(match(grid[selectionGrid[[i]],1],Mestimates$eval.points[[1]]),
                                      match(grid[selectionGrid[[i]],2],Mestimates$eval.points[[2]]))]+1E-10
      points=matrix(ncol=2,grid[selectionGrid[[i]],])
      if(length(selectionGrid[[i]])==0){points <- matrix(ncol=2,shapefile@polygons[[i]]@labpt)
      probs=1}
      if(npoints[i]>0){
        new=rbind(new,points[sample(1:nrow(points),size=max(0,npoints[i],na.rm=T),replace=T,prob=probs),])
      }
    }
    #recompute H
    if(adaptive==FALSE){
      H <- ks::Hpi(x = new, binned=TRUE)
    }
    if(adaptive==TRUE){
      H <- ks::Hpi(x = new, binned=TRUE)
      H <- sqrt(sqrt(H[1,1]*H[2,2]))
    }
    #recompute density
    if(adaptive==FALSE){
      Mestimates <- ks::kde(x=new, H=H,gridsize=c(length(gridx),length(gridy)),bgridsize=c(length(gridx),length(gridy)),
                            xmin=c(min(gridx),min(gridy)),xmax=c(max(gridx),max(gridy)),binned=TRUE)
    }
    if(adaptive==TRUE){
      counts <- plyr::count(new)
      MestimatesAd <- sparr::bivariate.density(data=counts[,c(1:2)],pilotH=H,res=length(gridx),xrange=range(gridx),
                                               yrange=range(gridy),adaptive=TRUE,
                                               comment=FALSE, counts=counts[,3])
      Mestimates$estimate=MestimatesAd$Zm
    }
    Mestimates$estimate[is.na(Mestimates$estimate)]=1E-96
    resultDensity[j,,] <- Mestimates$estimate
    resultX[j,,] <- new
    print(paste("Iteration:",j,"of", burnin+samples))
  }
  Mestimates$estimate=apply(resultDensity[-c(1:burnin),,],c(2,3),mean)
  est<-list(Mestimates=Mestimates,resultDensity=resultDensity,resultX=resultX,
            data=data, gridx=gridx, gridy=gridy, shapefile=shapefile,
            burnin=burnin, samples=samples, adaptive=adaptive)
  class(est) <- "bivshape"
  return(est)
}