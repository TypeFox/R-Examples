"simm.mba" <- function(date=1:100, x0=c(0,0), mu=c(0,0), sigma=diag(2),
                       id="A1", burst=id)
  {
      if (!inherits(date, "POSIXct")) {
          class(date) <- c("POSIXct", "POSIXt")
          attr(date, "tzone") <- ""
      }

      if (!is.matrix(sigma))
      stop("matrix expected")
    if (!all(dim(sigma)==2))
      stop("matrix 2*2 expected")
    if (sum((t(sigma)-sigma)^2)> 1e-7)
      stop("symetric matrix expected")
    if (any(eigen(sigma, symmetric=TRUE)$value< -1e-7))
      stop("positive matrix expected")

    n <- length(date)
    dt <- c(diff(unclass(date)),NA)

    dx <- c(rnorm(n-1,0,sqrt(dt[-n])),NA)
    dy <- c(rnorm(n-1,0,sqrt(dt[-n])),NA)
    W <- cbind(dx,dy)%*%sigma
    W[,1] <- W[,1] + mu[1]*dt
    W[,2] <- W[,2] + mu[2]*dt

    x <- c(x0[1], x0[1]+cumsum(W[-n,1]))
    y <- c(x0[2], x0[2]+cumsum(W[-n,2]))
    res <- as.ltraj(data.frame(x,y),date, id, burst, typeII=TRUE)
    return(res)
  }

