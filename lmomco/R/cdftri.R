"cdftri" <-
function(x,para) {
    if(! are.partri.valid(para)) return()
    MIN  <- para$para[1]
    MODE <- para$para[2]
    MAX  <- para$para[3]
    A <- (MAX  - MIN)
    B <- (MODE - MIN)
    C <- (MAX  - MODE)
    AB <- A*B
    AC <- A*C

    f <- vector(mode = "numeric", length=length(x))
    f[x >  MODE] <- 1 - (MAX-x[x > MODE]     )^2/(AC)
    f[x <  MODE] <-     (    x[x < MODE]-MIN)^2/(AB)
    f[x == MODE] <- B/A

    names(f) <- NULL
    f[! is.finite(f)] <- NA
    f[x < MIN] <- 0
    f[x > MAX] <- 1
    return(f)
}

