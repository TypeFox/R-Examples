# randompairs 15_06_12

randompairs <- function(nstim=10, nint=100, reverse=0.1, skew=FALSE) {
  if (nstim <= 26) {
    IDs <- sort(sample(letters, nstim))
  }
  if (nstim > 26 & nstim <= 325) {
    com <- combn(26, 2)[, -177]
    samplecom <- com[, sample(1:ncol(com), nstim)]
    IDs <- apply(samplecom, 2, function(x) letters[x])
    IDs <- sort(apply(IDs, 2, function(x) paste(x[1], x[2], sep = "")))
  }
  if (nstim > 325 & nstim <= 2601) {
    com <- combn(26, 3)
    samplecom <- com[, sample(1:ncol(com), nstim)]
    IDs <- apply(samplecom, 2, function(x) letters[x])
    IDs <- sort(apply(IDs, 2, function(x) paste(x[1], x[2], x[3], sep = "")))
  }


  if(!skew)  {
    xdata <- cbind(sample(IDs, nint*2, TRUE), sample(IDs, nint*2, TRUE))
  }
  if(skew)  {
    weights <- rnbinom(nstim, 10, 0.5)
    xdata <- cbind(sample(IDs, nint*2, TRUE, weights), sample(IDs, nint*2, TRUE, weights))
  }

  xdata <- xdata[ xdata[,1]!=xdata[,2], ]
  xdata <- xdata[1:nint, ]
  xdata <- t(apply(xdata, 1, sort))
  s <- sample(1:nrow(xdata), ceiling(nrow(xdata)*reverse))
  xdata[s, 1:2] <- xdata[s, 2:1]

  xdata <- data.frame(index=1:nrow(xdata), winner=xdata[,1], loser=xdata[,2])
  return(xdata)

}
