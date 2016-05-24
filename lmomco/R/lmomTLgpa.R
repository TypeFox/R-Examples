"lmomTLgpa" <-
function(para) {
    L <- seq(1,4)
    R <- seq(1,4)

    if(! are.pargpa.valid(para)) return()
    attributes(para$para) <- NULL

    XI <- para$para[1] 
    A  <- para$para[2]
    K  <- para$para[3]
    L[1] <- XI + A*(K+5)/((K+2)*(K+3))
    L[2] <- 6*A/((K+2)*(K+3)*(K+4))
    R[1] <- 0
    R[2] <- L[2]/L[1]
    R[3] <- 10*(1-K)/(9*(K+5))
    L[3] <- R[3]*L[2]
    R[4] <- 5*(K-1)*(K-2)/(4*(K+6)*(K+5))
    L[4] <- L[2]*R[4]
 
    z <- list(lambdas = L, ratios = R,
              trim = 1, leftrim = 1, rightrim=1,
              source = 'lmomTLgpa')
    return(z)
}
