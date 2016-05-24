ss.filter <- function(x,seasonal = FALSE, period = NULL,plot = TRUE,...)
{
  if (NCOL(x) > 1)
    stop("'x' must be a numeric vector or an univariate time series")
  if (any(!is.finite(x)))
    stop("missing values exist")
  n <- length(x)
  tn <- (1:n)/n
  xt <- x
  if (seasonal) {
    if (is.null(period))
      stop("'period' must be provided")
    if (period%%1 != 0 || period < 0)
      stop("'period' must be a positive integer")
    s.index <- rep_len(1:period,length.out = n)
    st <- numeric(n)
    for (i in 1:period) {
      st[s.index == i] <- mean(xt[s.index == i])
    }
    st <- st - mean(st)
    xt <- xt - st
  }
  tr.s <- smooth.spline(tn,xt,...)
  mt <- tr.s$y
  rt <- xt - mt
  result <- as.matrix(cbind(x, mt, if (seasonal) st, rt))
  colnames(result) <- c("data","trend",if (seasonal) "season","residual")
  if (plot) {
    if (seasonal) {
      op <- par(mfrow = c(2,2))
      plot(result[,1],type = "l",xlab = "time",ylab = "data",
           main = "original v.s fitted data")
      lines(result[,2] + result[,3], col = 2)
      plot(result[,2],type = "l", col = 3,xlab = "time",ylab = "estimated trend",
           main = "plot of estimated trend")
      plot(result[,3],type = "l", col = 4,xlab = "time",ylab = "seasonal index",
           main = "plot of seasonal index")
      plot(result[,4],type = "l", xlab = "time",ylab = "residuals",
           main = "plot of residuals")
      abline(h = 0, lty = 2, col = 5)
      par(op)
    }
    else {
      op <- par(mfrow = c(2,1))
      plot(result[,1],type = "l",xlab = "time",ylab = "data",
           main = "original v.s fitted data")
      lines(result[,2],col = 2)
      plot(result[,3],type = "l",xlab = "time",ylab = "residuals",
           main = "plot of residuals")
      abline(h = 0, lty = 2, col = 3)
      par(op)
    }
  }
  return(result) 
}