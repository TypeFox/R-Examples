"lmomkur" <-
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
              source = "lmomkur"
             )
    if(! are.parkur.valid(para)) return(NULL)
    attributes(para$para) <- NULL
    A <- para$para[1]
    B <- para$para[2]
    
    ao1 <- 1 + 1/A
    B1 <- beta(ao1,   B)
    B2 <- beta(ao1, 2*B)
    B3 <- beta(ao1, 3*B)
    B4 <- beta(ao1, 4*B)
    B5 <- beta(ao1, 5*B)
    
    tmp <- B1 - 2*B2
    z$L1   <- B*B1
    z$L2   <- B*tmp
    z$TAU3 <- (B1 -  6*B2 +  6*B3)/tmp
    z$TAU4 <- (B1 - 12*B2 + 30*B3 -  20*B4)/tmp
    z$TAU5 <- (B1 - 20*B2 + 90*B3 - 140*B4 + 70*B5)/tmp
    z$LCV  <- z$L2/z$L1
    z$L3 <- z$TAU3*z$L2
    z$L4 <- z$TAU4*z$L2
    z$L5 <- z$TAU5*z$L2
    z <- lmorph(z)
    return(z)
}
