"lmomlmrq" <-
function(para) {
    z <- list(lambdas=NA,
              ratios=NA,
              trim=1,
              leftrim=NULL,
              rightrim=NULL,
              source = "lmomlmrq")

    if(! are.parlmrq.valid(para)) return()
    attributes(para$para) <- NULL

    U <- para$para[1]
    A <- para$para[2]

    z$lambdas <- c(U, (A+3*U)/6, (A+U)/6, (A+U)/12, (A+U)/20, (A+U)/30)
    names(z$lambdas) <- NULL
    L2 <- z$lambdas[2] 
    z$ratios <- c(NA, L2/z$lambdas[1])
    names(z$ratios) <- NULL
    n <- length(z$lambdas)
    for(r in 3:n) {
       z$ratios[r] <- z$lambdas[r]/L2 
    }
    return(z)
}

