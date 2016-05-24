plot.transformedTrain <- function(x,
                                  which=1:5,
                                  main,
                                  caption=c("Uniform on Trans. Obs. Time Test",
                                    "Berman's Test",
                                    "Log Survivor Function",
                                    "Lag 1 Transformed Intervals",
                                    "Variance vs Mean",
                                    "Martingale vs Trans. Time"),
                                  ask=TRUE,
                                  ...
                                  ) {

  ## Check that x is a transformedTrain
  if (!inherits(x, "transformedTrain"))
    stop("use only with \"transformedTrain\" objects")

  Y <- seq(x)
  nbSpikes <- length(x)
  
  show <- logical(6)
  show[which] <- TRUE

  if (ask) {
    op <- par(ask = TRUE)
    on.exit(par(op))
  } else {
    if (sum(show)==2) layout(matrix(1:2,nrow=2))
    if (2<sum(show) && sum(show)<5) layout(matrix(1:4,nrow=2))
    if (sum(show)>=5) layout(matrix(1:6,nrow=3))
  }

  mainGiven <- !missing(main)
  if (show[1]) {
    slopeKS <- length(x)/max(x)
    plot(as.numeric(x),Y,type="n",
         xlab="Transformed time",
         ylab="Cumulative number of events",
         main=ifelse(mainGiven,main,caption[1]),
         sub=ifelse(mainGiven,caption[1],"")
         )
    abline(a=0,b=slopeKS)
    abline(a=1.36*sqrt(nbSpikes),slopeKS,lty=2)
    abline(a=-1.36*sqrt(nbSpikes),slopeKS,lty=2)
    abline(a=1.63*sqrt(nbSpikes),slopeKS,lty=2)
    abline(a=-1.63*sqrt(nbSpikes),slopeKS,lty=2)
    lines(as.numeric(x),Y,col=2,lwd=2)
  }

  lambda <- 1-exp(-diff(x))
  if (show[2]) {
    plot(c(0,1),c(0,1),type="n",
         xlab=expression(U[(k)]),
         ylab="Cumulative distribution",
         main=ifelse(mainGiven,main,caption[2]),
         sub=ifelse(mainGiven,caption[2],"")
         )
    abline(a=0,b=1)
    abline(a=1.36/sqrt(nbSpikes-1),1,lty=2)
    abline(a=-1.36/sqrt(nbSpikes-1),1,lty=2)
    abline(a=1.63/sqrt(nbSpikes-1),1,lty=2)
    abline(a=-1.63/sqrt(nbSpikes-1),1,lty=2)
    lines(sort(lambda),(1:(nbSpikes-1))/(nbSpikes-1),col=2,lwd=2)
  }

  if (show[3]) {
    isi <- diff(x)
    nI <- length(isi)
    Y <- (nI:1)/nI
    X <- sort(isi)
    Yth <- exp(-X)
    Y95p <- qbinom(0.975,nI,Yth)/nI
    Y95m <- qbinom(0.025,nI,Yth)/nI
    Y99p <- qbinom(0.995,nI,Yth)/nI
    Y99m <- qbinom(0.005,nI,Yth)/nI
    maxId <- max(which(Yth>0.001))
    plot(c(0,X[maxId]),c(0.001,1),type="n",
         xlab=expression(Y[(k)]),
         ylab="Survivor Fct",
         main=ifelse(mainGiven,main,caption[3]),
         sub=ifelse(mainGiven,caption[3],""),
         log="y"
         )
    lines(X,Y95p,lty=2)
    lines(X,Y95m,lty=2)
    lines(X,Y99p,lty=2)
    lines(X,Y99m,lty=2)
    lines(X,Y,col=2,lwd=2)
  }

  if (show[4]) {
    plot(lambda[-length(lambda)],
         lambda[-1],
         xlab=expression(U[k]),
         ylab=expression(U[k+1]),
         pch=3,
         main=ifelse(mainGiven,main,caption[4]),
         sub=ifelse(mainGiven,caption[4],"")
         )
  }

  if (show[5]) {
    plot(varianceTime(x),
         style="Ogata",
         xlab="Length of Trans. Time Window",
         main=ifelse(mainGiven,main,caption[5]),
         sub=ifelse(mainGiven,caption[5],"")
         )
  }

  if (show[6]) {
    X <- c(0,as.numeric(x))
    M <- 0:(length(X)-1) - X
    mTh <- -1.63*sqrt(length(X)-1)
    MTh <- 1.63*sqrt(length(X)-1)
    plot(X,
         M,
         type="n",
         xlab="Transformed time",
         ylab="Martingale",
         ylim=c(min(c(M,mTh)),max(c(M,MTh))),
         main=ifelse(mainGiven,main,caption[6]),
         sub=ifelse(mainGiven,caption[6],"")
         )
    abline(h=0)
    abline(h=1.36*sqrt(length(X)-1),lty=2)
    abline(h=-1.36*sqrt(length(X)-1),lty=2)
    abline(h=1.63*sqrt(length(X)-1),lty=2)
    abline(h=-1.63*sqrt(length(X)-1),lty=2)
    lines(X,M,col=2,lwd=2)
  }
}
