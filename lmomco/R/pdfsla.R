"pdfsla" <-
function(x,para) {
    if(! are.parsla.valid(para)) return()
    EPS   <- sqrt(.Machine$double.eps)
    SMALL <- 1E-6 # hacked into limit based on behavior of dnorm()
    U <- para$para[1]
    A <- para$para[2]

    Y <- (x - U)/A 
    DEL <- dnorm(0) - dnorm(Y)
    f <- DEL/Y^2
    names(f) <- NULL
    f[abs(DEL) <= EPS |     Y == 0    ] <- 1/(2*sqrt(2))/sqrt(pi)
    f[    DEL  == 0   & abs(Y) < SMALL] <- 1/(2*sqrt(2))/sqrt(pi)
    f[is.na(f)] <- 0 # decision Dec. 2015
    return(f)
}

