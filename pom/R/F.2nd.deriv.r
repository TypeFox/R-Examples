F.2nd.deriv <- function(pt, FUN, ...){

FUN <- match.fun( FUN )

d <- length(pt)   # number of dimensions
hess <- matrix( 0, nrow=d, ncol=d )

eps <- 10e-7
h <- (eps^(0.25))*pt


for( i in 1:d ){
    ei <- rep(0,d)
    ei[i] <- 1

    # compute diagonal element
    hess[i,i] <- (-FUN(pt+2*h*ei, ...) + 16*FUN(pt+h*ei, ...) - 30*FUN(pt, ...) + 16*FUN(pt-h*ei, ...) - FUN(pt-2*h*ei, ...)) / (12*h[i]*h[i])

    if( (i+1) <= d ){
        for( j in (i+1):d ){
            ej <- rep(0,d)
            ej[j] <- 1

            # compute off diagonal element
            hess[i,j] <- (FUN(pt+h*ei+h*ej, ...) - FUN(pt+h*ei-h*ej, ...) - FUN(pt-h*ei+h*ej, ...) + FUN(pt-h*ei-h*ej, ...)) / (4*h[i]*h[j])

            # Assume symetric
            hess[j,i] <- hess[i,j]

        }
    }
}


return(hess)

}

