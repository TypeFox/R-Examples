"lmomwei" <-
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
              source = "lmomwei"
             )
    if(! are.parwei.valid(para)) return()
    attributes(para$para) <- NULL

    ZETA <- para$para[1]
    B    <- para$para[2]
    D    <- para$para[3]

    K  <- 1/D
    A  <- B/D
    XI <- ZETA - B 
    gev.para <- list(type = 'gev', para = c(XI,A,K))

    z <- lmomgev(gev.para)
    z <- lmorph(z)
    z$L1   <- -z$L1
    z$LCV  <- -z$LCV
    z$TAU3 <- -z$TAU3
    z$TAU5 <- -z$TAU5
    z$L3   <- -z$L3
    z$L5   <- -z$L5

    z$source <- "lmomwei"
    z <- lmorph(z)
    return(z)
}

