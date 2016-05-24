"quatri" <-
function(f,para, paracheck = TRUE) {
    if(! check.fs(f)) return()
    if(paracheck == TRUE) {
       if(! are.partri.valid(para)) return()
    }

    MIN  <- para$para[1]
    MODE <- para$para[2]
    MAX  <- para$para[3]
    A <- (MAX  - MIN)
    B <- (MODE - MIN)
    C <- (MAX  - MODE)
    AB <- A*B
    AC <- A*C

    P <- B/A # The F value at the mode

    x <- vector(mode = "numeric", length=length(f))
    x[f >  P] <- MAX - sqrt((1-f[f > P])*AC)
    x[f <  P] <- MIN + sqrt((  f[f < P])*AB)
    x[f == P] <- MODE

    return(x)
}

