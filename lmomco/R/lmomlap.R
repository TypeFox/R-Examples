"lmomlap" <-
function(para) {

    z <- list(lambdas = c(NA),
              ratios  = c(0),
              trim = 1,
              leftrim = NULL,
              rightrim = NULL,
              source = "lmomlap"
             )

    if(! are.parlap.valid(para)) return()
    attributes(para$para) <- NULL

    XI <- para$para[1]
    A  <- para$para[2]
    z$lambdas[1] <- XI
    z$lambdas[2] <- 3*A/4
    z$lambdas[3] <- 0
    z$lambdas[4] <- 17*A/96
    z$lambdas[5] <- 0
    z$lambdas[6] <- 31/480*A; # 0.06458333*A

    z$ratios[2]  <- z$lambdas[2]/z$lambdas[1]
    for(i in 3:6) {
       z$ratios[i] <- z$lambdas[i]/z$lambdas[2]
    }
    return(z)
}


