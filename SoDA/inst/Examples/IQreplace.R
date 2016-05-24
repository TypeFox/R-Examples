newIQ <- function(nData = 1000, probs = seq(0, 1, 0.25)) 
  list(nData = nData, probs = probs,
       dataBuf = numeric(0), qBuf = numeric(0))

`addData<-` <- function(IQ, update = FALSE, value) {
    n <- length(value);
    if(update || (n + length(IQ$dataBuf) > IQ$nData))
      recompute(IQ, value)
    else {
        IQ$dataBuf <- c(IQ$dataBuf, value)
        IQ
    }
}

recompute  <- function(IQ, newdata = numeric(0)) {
    IQ$qBuf <- doQuantile(qBuf, c(IQ$dataBuf, newdata), IQ$probs)
    IQ$dataBuf <- numeric(0)
    IQ
}

getQ <- function(IQ) {
    if(length(IQ$dataBuf) > 0)
      doQuantile(IQ$qbuf, IQ$dataBuf, IQ$probs)
    else
      IQ$qBuf
}
