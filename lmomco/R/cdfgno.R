"cdfgno" <-
function(x,para) {
    # Error function from R documentation
    erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1
    RTHALF <- 0.707106781186547524 # cos(pi/4)

    if(! are.pargno.valid(para)) return()
    XI <- para$para[1]
    A  <- para$para[2]
    K  <- para$para[3]

    Y <- (x-XI)/A
    if(K != 0) {
       ARG <- 1-K*Y
       ops <- options(warn = -1)
       Y <- -log(ARG)/K
       ops <- options(ops)
    }
    f <- 0.5+0.5*erf(Y*RTHALF)
    if(K < 0) {
       f[!is.finite(f)] <- 0
    } else if(K > 0) {
       f[!is.finite(f)] <- 1
    }
    names(f) <- NULL
    return(f)
}

