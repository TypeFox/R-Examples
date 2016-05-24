"partri" <-
function(lmom,checklmom=TRUE) {

    if(length(lmom$L1) == 1) { # convert to named L-moments
      lmom <- lmorph(lmom)     # nondestructive conversion!
    }
    if(checklmom & ! are.lmom.valid(lmom)) {
      warning("L-moments are invalid")
      return()
    }
    
    para <- c(NA, NA, NA)
    names(para) <- c("min","mode","max")

    S  <- 30/7
    L1 <- lmom$lambdas[1]
    L2 <- lmom$lambdas[2]
    T3 <- lmom$ratios[3]
    Xg <- L1 + S*L2
    Ng <- L1 - S*L2

    EPS <- 1E-7
    maxTau3 <- 0.14285710
    #print(T3)
    #print(T3+EPS)
    if(T3+EPS < -maxTau3 | T3-EPS > maxTau3) {
       warning("abs(L-skew) > ",maxTau3," and is numerically incompatible with the distribution")
       return(list(type = 'tri', para=para, obj.val="abs(Tau3) too big", source="partri"))
    } else if(abs(T3) >= maxTau3) {
       T3 <- sign(T3)*maxTau3 # if T3 is 0.1 ppm close to the maximum, set it as such
    }

    Og <- (Xg + Ng)/2

    init.par <- c(Ng, Og, Xg)
    #print(init.par)

    "afunc" <- function(par) {
       z <- list(para=par, type="tri")
       lmr <- lmomtri(z, paracheck=FALSE)
       err <- ((lmr$lambdas[1] - L1)/L1)^2 +
              ((lmr$lambdas[2] - L2)/L2)^2 +
              ((lmr$ratios[3]  - T3)   )^2
       return(err)
    }

    opt <- NULL
    try(opt <- optim(init.par, afunc))
    if(is.null(opt)) {
       warning("Could not simultaneously optimize the three parameters, returning NULL")
       return(NULL)
    }
    if(opt$convergence != 0) {
       warning("Convergence problems, printing out results")
       print(opt)
    }
    if(opt$par[2] > opt$par[3] | opt$par[2] < opt$par[1]) {
       message("Likely just numerical problem encountered if the distribution is at or ",
               "nearly at right triangle but the mode is incompatible with one or the ",
               "limits, proceeding to sort the parameter estimates")
       message("Numerical estimate: ",paste(opt$par, collapse=" "))
       para <- sort(opt$par)
       message("               New: ",paste(para,    collapse=" "))
    } else {
       para <- opt$par
    }
    if(T3 >  maxTau3) { # snap to a right triangle
       min <- (para[1]+para[2])/2
       para[1] <- para[2] <- min
    }
    if(T3 < -maxTau3) { # snap to a right triangle
       max <- (para[2]+para[3])/2
       para[2] <- para[3] <- max
    }
    names(para) <- c("min","mode","max")

    return(list(type = 'tri', para = para, obj.val=opt$value, source="partri"))
}

