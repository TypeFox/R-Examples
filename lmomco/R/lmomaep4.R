"lmomaep4" <-
function(para, paracheck=TRUE, t3t4only=FALSE) {

    if(paracheck && ! are.paraep4.valid(para)) {
       warning("Parameters are invalid")
       return()
    }
    attributes(para$para) <- NULL



    z <- list(lambdas  = rep(NA,4),
              ratios   = rep(NA,4),
              trim     = 0,
              leftrim  = 0,
              rightrim = 0,
              source = "lmomaep4"
             )


    U <- para$para[1]
    A <- para$para[2]
    K <- para$para[3]
    H <- para$para[4]

    # Delicado and Goria (2008, p.6)
    KmK <- 1/K - K
    H2H1 <- exp(lgamma(2/H) - lgamma(1/H))
    Ihalf <- pbeta(1/2, shape1=1/H, shape2=2/H)

    KK    <- K*K
    KKK   <- KK*K
    KKKK  <- KKK*K
    KKKKK <- KKKK*K

    L1 <- U + A * KmK * H2H1
    z$lambdas[1] <- L1

    L2a <- -K * KmK^2              / (1+KK)
    L2b <-  2 * KK * (1/KKK + KKK) / (1+KK)^2 * Ihalf
    L2  <-  (L2a + L2b) * H2H1
    z$lambdas[2] <- L2 * A
    z$ratios[2]  <- z$lambdas[2]/L1

    "delf" <- function(t) {
       delfa <- t^(1/H - 1)
       delfb <- (1 - t)^(2/H - 1)
       It <- pbeta((1-t)/(2-t), shape1=1/H, shape2=3/H)
       return(delfa*delfb*It)
    }
    DELTAint <- NULL
    try( DELTAint <- integrate(delf, 0, 1/2) )
    if(is.null(DELTAint)) return(z)
    DELTA <- DELTAint$value
    DELTA <- DELTA/beta(1/H, 2/H)

    L3a <-  1 * KmK * (KKKK - 4*KK + 1)   / (1+KK)^2
    L3b <- -6 * KKK * KmK * (1/KKK + KKK) / (1+KK)^3 * Ihalf
    L3c <-  6 * (1+KKKK) * KmK            / (1+KK)^2 * DELTA
    L3  <- (L3a + L3b + L3c) * H2H1
    z$lambdas[3] <- L3 * A
    z$ratios[3]  <- L3/L2

    "delf1y" <- function(y) {
         return(sapply(y,
            function(ay) {
               "delf1z" <- function(z) {
                  t1 <- ay^(1/H - 1) * (1 - ay)^(2/H - 1)
                  t2 <-  z^(1/H - 1) *  (1 - z)^(3/H - 1)
                  j <- (1-z)*(1-ay)
                  j <- j / (1+j)
                  It <- pbeta(j, shape1=1/H, shape2=4/H)
                  return(t1*t2*It)
              }
              try( DELTintz <- integrate(delf1z, 0, (1-ay)/(2-ay) ) )
              return(DELTintz$value)
            }))
    }
    DELTAint1 <- NULL
    try( DELTAint1 <- integrate(delf1y, 0, 1/2) )
    if(is.null(DELTAint1)) return(z)
    DELTA1 <- DELTAint1$value
    DELTA1 <- DELTA1 / (beta(1/H, 2/H) * beta(1/H, 3/H))

    L4a <-  -1 * K * KmK^2 * (KKKK - 8*KK + 1)          / (1+KK)^3
    L4b <-  12 * KK * (KKK + 1/KKK) * (KKKK - 3*KK + 1) / (1+KK)^4 * Ihalf
    L4c <- -30 * KKK * KmK^2 * (1/KK + KK)              / (1+KK)^3 * DELTA
    L4d <-  20 * KKKK * (1/KKKKK + KKKKK)               / (1+KK)^4 * DELTA1
    L4  <- (L4a + L4b + L4c + L4d) * H2H1
    z$lambdas[4] <- L4 * A
    z$ratios[4]  <- L4/L2

    if(t3t4only) {
      return(list(T3=z$ratios[3], T4=z$ratios[4]))
    } else {
      return(z)
    }
}
