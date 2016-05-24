"simm.crw" <- function(date=1:100, h = 1, r = 0,
                       x0=c(0,0), id="A1", burst=id,
                       typeII=TRUE)
{
    if (typeII) {
        if (!inherits(date, "POSIXct")) {
            class(date) <- c("POSIXct", "POSIXt")
            attr(date, "tzone") <- ""
        }
    }
    n <- length(date)
    dt <- c(diff(unclass(date)))
    if (all(dt-dt[1]>1e-7))
        stop("the time lag between relocations should be constant")

    ang<-rwrpnorm(n-2,0,r)
    if (h>0) {
        v=sqrt(dt)*rchi(n-1) * h
    } else {
        v=-h
    }
    ang=cumsum(c(runif(1,0,2*pi),ang))
    si=c(x0[2], x0[2]+cumsum(v*sin(ang)))
    co=c(x0[1], x0[1]+cumsum(v*cos(ang)))
    res <- as.ltraj(data.frame(co,si),date, id, burst,
                    typeII=typeII)
    return(res)
}

