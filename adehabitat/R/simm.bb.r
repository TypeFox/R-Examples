"simm.bb" <- function(date=1:100, begin = c(0,0),
                      end=begin, id="A1", burst=id)
  {
    class(date) <- c("POSIX","POSIXct")
    n <- length(date)
    dt <- c(diff(unclass(date)),NA)

    dx <- c(rnorm(n-1,0,sqrt(dt[-n])),NA)
    dy <- c(rnorm(n-1,0,sqrt(dt[-n])),NA)
    W <- cbind(dx,dy)

    xtmp <- begin
    x <- begin
    for (i in 2:(n-1)) {
      a <- diag(1/(date[n]-date[i-1]),2)
      dX <- as.numeric((a%*%(end - xtmp)*dt[i-1] + W[i-1,]))
      x <- rbind(x,xtmp+dX)
      xtmp <- xtmp+dX
    }
    x <- rbind(x,end)

    res <- as.ltraj(data.frame(x[,1],x[,2]),date, id, burst, typeII=TRUE)
    return(res)

  }

