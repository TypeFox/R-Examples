"lmomgum" <-
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
              source = "lmomgum"
             )
    #  ARRAY ZMOM CONTAINS THE L-MOMENT RATIOS OF THE STANDARD 
    #  GUMBEL DISTRIBUTION (XI=0, ALPHA=1). 
    #  ZMOM(1) IS EULER'S CONSTANT, ZMOM(2) IS LOG(2). 
    #
    ZMOM <- c(0.577215664901532861,
              0.693147180559945309, 
              0.169925001442312363,
              0.150374992788438185, 
              0.558683500577583138e-1)

    if(! are.pargum.valid(para)) return()
    attributes(para$para) <- NULL

    A <- para$para[2] 
    z$L1 <- para$para[1] + A*ZMOM[1] 
    z$L2 <- A*ZMOM[2] 
    z$TAU3 <- ZMOM[3]
    z$TAU4 <- ZMOM[4]
    z$TAU5 <- ZMOM[5]
    z$LCV  <- z$L2/z$L1
    z$L3   <- z$TAU3*z$L2
    z$L4   <- z$TAU4*z$L2
    z$L5   <- z$TAU5*z$L2
    z <- lmorph(z)
    return(z)
}

