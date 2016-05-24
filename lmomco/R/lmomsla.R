"lmomsla" <-
function(para) {

    if(! are.parsla.valid(para)) return()
    attributes(para$para) <- NULL
    
    L1 <- para$para[0]
    L2 <- 0.9368627*para$para[2]
    ratios  <- c(NA, L2/L1, 0,    0.3042045, 0,    0.1890072)
    lambdas <- c(L1,    L2, 0, ratios[4]*L2, 0, ratios[6]*L2)
    
    z <- list(lambdas=lambdas, ratios=ratios,
              trim=1, leftrim=1, rightrim=1, source='lmomsla')
    return(z)
}

