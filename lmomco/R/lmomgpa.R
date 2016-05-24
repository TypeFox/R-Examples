"lmomgpa" <-
function(para) {
    z <- list(L1   = NULL,
              L2   = NULL,
              TAU3 = NULL,
              TAU4 = NULL,
              TAU5 = NULL,
              LCV  = NULL,
              L3   = NULL,
              L4   = NULL,
              L5   = NULL,
              source = "lmomgpa"
             )

    if(! are.pargpa.valid(para)) return()
    attributes(para$para) <- NULL

    XI <- para$para[1] 
    A  <- para$para[2] 
    K  <- para$para[3] 

    #  LAMBDA-1
    Y <- 1/(1+K)
    z$L1 <- XI+A*Y

    #  LAMBDA-2
    Y <- Y/(2+K)
    z$L2 <- A*Y

    #  HIGHER MOMENTS
    x <- vector(mode="numeric",length=5)
    Y <- 1
    for(m in seq(3,5)) {
      AM   <- m-2
      Y    <- Y*(AM-K)/(m+K)
      x[m] <- Y
    }
    z$TAU3 <- x[3]
    z$TAU4 <- x[4]
    z$TAU5 <- x[5]
    z$LCV  <- z$L2/z$L1
    z$L3   <- z$TAU3*z$L2
    z$L4   <- z$TAU4*z$L2
    z$L5   <- z$TAU5*z$L2
    z <- lmorph(z)
    return(z)
}

