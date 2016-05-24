"ht2" <-
function (data) 
{
    a <-2
    n <- length(data)
    nhalf <- n/2
    J <- logb(n, 2)
    res <- data
res1<-res2<-NULL
    sm <- rep(0, nhalf)
    det <- sm
    for (i in 1:J) {
        sm[1:nhalf] <- (res[2 * (1:nhalf) - 1] + res[2 * (1:nhalf)])/a
        det[1:nhalf] <- (res[2 * (1:nhalf) - 1] - res[2 * (1:nhalf)])/a
res1<-c(res1,det[1:nhalf])
res2<-c(res2,sm[1:nhalf])
        res[1:nhalf] <- sm[1:nhalf]
        res[(nhalf + 1):n] <- det[1:nhalf]
        n <- n/2
        nhalf <- nhalf/2
        sm <- 0
        det <- 0
    }

    return(list(s=res2,d=res1))

#res
}

