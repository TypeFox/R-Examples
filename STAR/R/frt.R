"%frt%" <- function(refTrain,testTrain) frt(refTrain,testTrain)

frt <- function(refTrain,
                testTrain
                ) {

  class(refTrain) <- NULL
  class(testTrain) <- NULL
  refTrain <- refTrain[!is.na(refTrain)]
  testTrain <- testTrain[!is.na(testTrain)]
  result <- sapply(refTrain,
                   function(t) {
                     goodOnes <- testTrain >= t
                     if (sum(goodOnes) == 0) return(NA)
                     min(testTrain[goodOnes])-t
                   }
                   )
  result <- result[!is.na(result)]
  class(result) <- "frt"
  result
}

plot.frt <- function(x,
                     which=1:2,
                     main,
                     caption=c("Log Survivor Function",
                       "Berman's Test"),
                     ask=TRUE,
                     ...) {

  show <- logical(2)
  show[which] <- TRUE

  if (ask) {
    op <- par(ask = TRUE)
    on.exit(par(op))
  } else {
    if (sum(show)==2) layout(matrix(1:2,nrow=2))
  }

  X <- sort(x)
  nI <- length(X)
  mainGiven <- !missing(main)
  
  if (show[1]) {
    Y <- (nI:1)/nI
    YTh <- exp(-X)
    p95 <- qbinom(0.975,nI,YTh)/nI
    m95 <- qbinom(0.025,nI,YTh)/nI
    p99 <- qbinom(0.995,nI,YTh)/nI
    m99 <- qbinom(0.005,nI,YTh)/nI
    plot(X,Y,
         type="n",log="y",
         xlab="Rescaled Forward Recurrence Time",
         ylab="Survivor Fct",
         main=ifelse(mainGiven,main,caption[1]),
         sub=ifelse(mainGiven,caption[1],""),
         ...)
    lines(X,YTh)
    lines(X,p95,lty=2)
    lines(X,m95,lty=2)
    lines(X,p99,lty=2)
    lines(X,m99,lty=2)
    lines(X,Y,col=2,lwd=2)
  }

  lambda <- 1-exp(-x)
  
  if (show[2]) {
    plot(c(0,1),c(0,1),type="n",
         xlab=expression(U[(k)]),
         ylab="Cumulative Distribution",
         main=ifelse(mainGiven,main,caption[2]),
         sub=ifelse(mainGiven,caption[2],""),
         ...
         )
    abline(a=0,b=1)
    abline(a=1.36/sqrt(nI),1,lty=2)
    abline(a=-1.36/sqrt(nI),1,lty=2)
    abline(a=1.63/sqrt(nI),1,lty=2)
    abline(a=-1.63/sqrt(nI),1,lty=2)
    lines(sort(lambda),(1:(nI))/nI,col=2,lwd=2)
  }

}

summary.frt <- function(object,...) {

  lambda <- sort(1-exp(-object))
  nI <- length(lambda)
  M <- (1:(nI))/nI - lambda
  in95 <- all(-1.36/sqrt(nI) <= M & M <= 1.36/sqrt(nI))
  in99 <- all(-1.63/sqrt(nI) <= M & M <= 1.63/sqrt(nI))
  c(in95=in95,in99=in99)
}
