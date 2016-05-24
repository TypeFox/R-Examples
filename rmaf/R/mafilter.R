ma.filter <- function (x, q = NULL, seasonal = FALSE, period = NULL,plot = TRUE)
{
  if (NCOL(x) > 1)
    stop("'x' must be a numeric vector or an univariate time series")
  if (any(!is.finite(x)))
    stop("missing values are not allowed")
  n <- length(x)
  if (n < 1L)
    stop("invalid length of 'x'") 
  trendest <- function(xt,q) {
    xt <- as.vector(xt)
    n <- as.integer(length(xt))
    mhat <- NULL
    if (is.null(q))
      q.n <- ifelse(qn(xt) < n, as.integer(min(qn(xt), n - qn(xt))), floor(n^(4/5)/2))
    else q.n <- q
    for (j in 1:n) {
      if (j >= q.n + 1 && j <= n - q.n) {
        mhat[j] <- sum(xt[(j - q.n):(q.n + j)])/(2 * q.n + 1)
      }
      else if (j < q.n + 1) {
        mhat[j] <- ((4 * q.n^2 - 4 * q.n * j + 6 * q.n + 4 * j^2 - 6 * j + 2) * 
                      sum(xt[1:(q.n + j)]) - 6 * (q.n - j + 1) * 
                      sum(c((-j + 1):q.n) * xt[1:(q.n + j)]))/((q.n + j) * ((q.n + j)^2 - 1))
      }
      else if (j > n - q.n) {
        mhat[j] <- ((4 * (n - j)^2 + 4 * q.n * (q.n + j - n) + 2 * (n + q.n - j)) * 
                      sum(xt[(j - q.n):n]) + 6 * (q.n - n + j) * 
                      sum(c(-q.n:(n - j)) * xt[(j - q.n):n]))/
                      ((n + q.n - j + 2) * (n + q.n - j + 1) * (n + q.n - j))
      }
    }
    return(mhat)
  }
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
  mt <- trendest(xt,q)
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