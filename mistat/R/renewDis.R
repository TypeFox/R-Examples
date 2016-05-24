renewDis <- function(ttf, ttr, time, n, printSummary=TRUE){
  ttf <- as.numeric(ttf)
  ttr <- as.numeric(ttr)
  time <- as.numeric(time)
  n <- as.integer(n)
  if(length(time) != 1 || length(n) != 1)
    stop("time and n should be of lenght 1")
  res <- as.matrix(rep(as.numeric(NA), n), ncol=1)
  for(i in 1:n){
    tt <- 0
    it <- 1
    while(tt < time){
      tt <- sum(c(tt, sample(ttf, 1)))
      tt <- sum(c(tt, sample(ttr, 1)))
      it <- it + 1
    }
    res[i,1] <- it
  }
  if(printSummary){
    cat(paste(" The estimated MEAN NUMBER Of RENEWALS is", round(mean(res[,1]), 2)))
    cat("\n")
    cat("number of renewals EBD\n")
    print(summary(c(res)))
  } 
  invisible(as.numeric(res))
}