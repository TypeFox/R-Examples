"simm.brown" <- function(date=1:100, x0=c(0,0), h = 1, id="A1", burst=id)
{
    if (!inherits(date, "POSIXct")) {
        class(date) <- c("POSIXct", "POSIXt")
        attr(date, "tzone") <- ""
    }
    y0 <- x0[2]
    x0 <- x0[1]
    n <- length(date)
    dt <- c(diff(unclass(date)),NA)
    dx <- c(rnorm(n-1,0,sqrt(dt[-n])*h),NA)
    dy <- c(rnorm(n-1,0,sqrt(dt[-n])*h),NA)
    x <- c(x0, x0+cumsum(dx[-n]))
    y <- c(y0, y0+cumsum(dy[-n]))
    res <- as.ltraj(data.frame(x,y),date, id, burst, typeII=TRUE)
    return(res)
}

