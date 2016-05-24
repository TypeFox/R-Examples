"lmomcau" <-
function(para) {
    L <- seq(1,4)
    R <- seq(1,4)
    attributes(para$para) <- NULL

    XI <- para$para[1]
    A  <- para$para[2]
    L[1] <- XI 
    L[2] <- 0.698*A 
    R[1] <- 0
    R[2] <- L[2]/L[1]
    R[3] <- 0
    L[3] <- 0
    R[4] <- 0.343
    L[4] <- L[2]*R[4]
 
    z <- list(lambdas = L, ratios = R, trim=1, leftrim=1, rightrim=1, source = "lmomcau")
    return(z)
}
