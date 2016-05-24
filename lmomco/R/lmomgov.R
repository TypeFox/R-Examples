"lmomgov" <-
function(para) {

    if(! are.pargov.valid(para)) {
       warning("Parameters are invalid")
       return()
    }

    names(para$para) <- NULL
    U <- para$para[1]
    A <- para$para[2]
    B <- para$para[3]

    L1 <- U + 2*A/(B+2)
    L2 <- 2*B*A/((B+2)*(B+3))
    T3 <- (B-2)/(B+4)
    T4 <- (B-5)*(B-1)/((B+4)*(B+5))

    T2 <- L2/L1; L3 <- L2*T3; L4 <- L2*T4
    z <- list(lambdas=c(L1,L2,L3,L4),
              ratios=c(NA, T2, T3, T4),
              source="lmomgov")
    return(z)
}


