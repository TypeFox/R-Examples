"simm.mou" <-
function(date=1:100, b = c(0,0), a = diag(0.5,2),
                     x0=b,  sigma=diag(2),
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

    xtmp <- x0
    x <- x0
    for (i in 2:n) {
      dX <- as.numeric((a%*%(b - xtmp)*dt[i-1] + W[i-1,]))
      x <- rbind(x,xtmp+dX)
      xtmp <- xtmp+dX
    }

    res <- as.ltraj(data.frame(x[,1],x[,2]),date, id, burst, typeII=TRUE)
    return(res)
  }

