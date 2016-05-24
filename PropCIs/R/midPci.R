midPci <-
function(x,n,conf.level){
    alpha <- 1 - conf.level
    pp<-seq(0.0001, 1 , 0.0005)
    uplim<-1
    lowlim<-0
    if (x == 0)
        uplim <- 1-alpha^(1/n)
    if (x == n)
        lowlim <- (alpha)^(1/n)
    if (x>0 & x<n){
        a2 <- 0.5*pbinom(x-1, n , pp,lower.tail = TRUE) +
            0.5*pbinom(x, n , pp, lower.tail = TRUE, log.p = FALSE)
        uplim=pp[ max(which(a2>(alpha/2))) ]
        lowlim=pp[ min(which(a2<(1-alpha/2))) ]
    }
    cint <- c(lowlim, uplim)
    attr(cint, "conf.level") <- conf.level
    rval <- list(conf.int = cint)
    class(rval) <- "htest"
    return(rval)
}

