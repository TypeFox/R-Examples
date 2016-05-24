BiCopPar2Tau <- function(family, par, par2 = 0, obj = NULL) {
    ## extract family and parameters if BiCop object is provided
    if (missing(family))
        family <- NA
    if (missing(par))
        par <- NA
    if (!is.null(obj)) {
        stopifnot(class(obj) == "BiCop")
        family <- obj$family
        par <- obj$par
        par2 <- obj$par2
    }
    if (class(family) == "BiCop") {
        # for short hand usage extract from family
        obj <- family
        family <- obj$family
        par <- obj$par
        par2 <- obj$par2
    }
    
    ## sanity checks for family and parameters
    if (is.na(family) || is.na(par))
        stop("Provide either 'family' and 'par' or 'obj'")
    if (length(family) != 1) {
        stop("Input for family has to be a scalar/integer.")
    }
    if (!(family %in% c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
                        13, 14, 16, 17, 18, 19, 20, 
                        23, 24, 26, 27, 28, 29, 30,
                        33, 34, 36, 37, 38, 39, 40, 
                        41, 42, 51, 52, 61, 62, 71, 72,
                        104, 114, 124, 134, 204, 214, 224, 234))) {
        stop("Copula family not implemented.")
    }
    
    if (missing(par)) {
        stop("'par' not set.")
    }
    
    if (length(par2) > 1 && length(par) != length(par2)) {
        stop("Input for 'par' and 'par2 has to be vectors of same length.")
    }
    
    if (family %in% c(7, 8, 9, 10,
                      17, 18, 19, 20, 
                      27, 28, 29, 30,
                      37, 38, 39, 40, 
                      42, 52, 62, 72,
                      104, 114, 124, 134,
                      204, 214, 224, 234) && par2 == 0) {
        stop("For BB1, BB6, BB7, BB8 and Tawn copulas, 'par2' must be set.")
    }
    
    if ((family == 1 || family == 2) && any(abs(par[1]) >= 1)) {
        stop("The (first) parameter of the Gaussian and t-copula has to be in the interval (-1,1).")
    }
    # if(family==2 && par2<=2) stop('The degrees of freedom parameter of the t-copula
    # has to be larger than 2.')
    if ((family == 3 || family == 13) && any(par <= 0)) {
        stop("The parameter of the (survival) Clayton copula has to be positive.")
    }
    if ((family == 4 || family == 14) && any(par < 1)) {
        stop("The parameter of the (survival) Gumbel copula has to be in the interval [1,oo).")
    }
    if ((family == 6 || family == 16) && any(par <= 1)) {
        stop("The parameter of the (survival) Joe copula has to be in the interval (1,oo).")
    }
    if (family == 5 && any(par == 0)) {
        stop("The parameter of the Frank copula has to be different from 0.")
    }
    if ((family == 7 || family == 17) && any(par <= 0)) {
        stop("The first parameter of the BB1 copula has to be positive.")
    }
    if ((family == 7 || family == 17) && any(par2 < 1)) {
        stop("The second parameter of the BB1 copula has to be in the interval [1,oo).")
    }
    if ((family == 8 || family == 18) && any(par <= 0)) {
        stop("The first parameter of the BB6 copula has to be in the interval [1,oo).")
    }
    if ((family == 8 || family == 18) && any(par2 < 1)) {
        stop("The second parameter of the BB6 copula has to be in the interval [1,oo).")
    }
    if ((family == 9 || family == 19) && any(par < 1)) {
        stop("The first parameter of the BB7 copula has to be in the interval [1,oo).")
    }
    if ((family == 9 || family == 19) && any(par2 <= 0)) {
        stop("The second parameter of the BB7 copula has to be positive.")
    }
    if ((family == 10 || family == 20) && any(par < 1)) {
        stop("The first parameter of the BB8 copula has to be in the interval [1,oo).")
    }
    if ((family == 10 || family == 20) && any(par2 <= 0 || par2 > 1)) {
        stop("The second parameter of the BB8 copula has to be in the interval (0,1].")
    }
    if ((family == 23 || family == 33) && any(par >= 0)) {
        stop("The parameter of the (90/270 degree) rotated Clayton copula has to be negative.")
    }
    if ((family == 24 || family == 34) && any(par > -1)) {
        stop("The parameter of the (90/270 degree) rotated Gumbel copula has to be in the interval (-oo,-1].")
    }
    if ((family == 26 || family == 36) && any(par >= -1)) {
        stop("The parameter of the (90/270 degree) rotated Joe copula has to be in the interval (-oo,-1).")
    }
    if ((family == 27 || family == 37) && any(par >= 0)) {
        stop("The first parameter of the (90/270 degree) rotated BB1 copula has to be negative.")
    }
    if ((family == 27 || family == 37) && any(par2 > -1)) {
        stop("The second parameter of the (90/270 degree) rotated BB1 copula has to be in the interval (-oo,-1].")
    }
    if ((family == 28 || family == 38) && any(par >= 0)) {
        stop("The first parameter of the rotated BB6 copula has to be in the interval (-oo,-1].")
    }
    if ((family == 28 || family == 38) && any(par2 > -1)) {
        stop("The second parameter of the (90/270 degree) rotated BB6 copula has to be in the interval (-oo,-1].")
    }
    if ((family == 29 || family == 39) && any(par > -1)) {
        stop("The first parameter of the (90/270 degree) rotated BB7 copula has to be in the interval (-oo,-1].")
    }
    if ((family == 29 || family == 39) && any(par2 >= 0)) {
        stop("The second parameter of the (90/270 degree) rotated BB7 copula has to be negative.")
    }
    if ((family == 30 || family == 40) && any(par > -1)) {
        stop("The first parameter of the (90/270 degree) rotated BB8 copula has to be in the interval (-oo,-1].")
    }
    if ((family == 30 || family == 40) && any(par2 >= 0 || par2 < (-1))) {
        stop("The second parameter of the (90/270 degree) rotated BB8 copula has to be in the interval [-1,0).")
    }
    if ((family == 41 || family == 51) && any(par <= 0)) {
        stop("The parameter of the (survival) reflection asymmetric copula has to be positive.")
    }
    if ((family == 61 || family == 71) && any(par >= 0)) {
        stop("The parameter of the (90/270 degree) rotated reflection asymmetric copula has to be negative.")
    }
    if (family == 42) {
        a <- par
        b <- par2
        limA <- (b - 3 - sqrt(9 + 6 * b - 3 * b^2))/2
        if (any(abs(b) > 1)) 
            stop("The second parameter of the two-parametric asymmetric copulas has to be in the interval [-1,1]")
        if (any(a > 1 || a < limA)) 
            stop("The first parameter of the two-parametric asymmetric copula has to be in the interval [limA(par2),1]")
    }
    if ((family == 104 || family == 114 || family == 204 || family == 214) && any(par < 1)) {
        stop("Please choose 'par' of the (180 degree rotated) Tawn copula in [1,oo).")
    }
    if ((family == 104 || family == 114 || family == 204 || family == 214) && any(par2 < 0 || par2 > 1)) {
        stop("Please choose 'par2' of the (180 degree rotated) Tawn copula in [0,1].")
    }
    if ((family == 124 || family == 134 || family == 224 || family == 234) && any(par > -1)) {
        stop("Please choose 'par' of the (90/270 degree) rotated Tawn copula in (-oo,-1].")
    }
    if ((family == 124 || family == 134 || family == 224 || family == 234) && any(par2 < 0 || par2 > 1)) {
        stop("Please choose 'par2' of the (90/270 degree) rotated Tawn copula in [0,1].")
    }
    
    ## calculation of tau(s) depending on pair-copula family
    if (family == 0) {
        tau <- rep(0, times = length(par))
    } else if (family == 1 | family == 2) {
        tau <- 2/pi * asin(par)
    } else if (family == 3 || family == 13) {
        tau <- par/(par + 2)
    } else if (family == 4 || family == 14) {
        tau <- 1 - 1/par
    } else if (family == 5) {
        f <- function(x) x/(exp(x) - 1)
        fu <- function(x) integrate(f, lower = 0, upper = x)$value
        fl <- function(x) integrate(f, lower = x, upper = 0)$value
        if (any(par > 0)) {
            tau <- 1 - 4/par + 4/par^2 * sapply(par, fu)
        } else {
            tau <- 1 - 4/par - 4/par^2 * sapply(par, fl)
        }
    } else if (family == 6 || family == 16) {
        # tau = 1 + 4/par^2 * integrate(function(x) log(x)*x*(1-x)^(2*(1-par)/par), 0,
        # 1)$value
        param1 <- 2/par + 1
        tem <- digamma(2) - digamma(param1)
        tau <- 1 + tem * 2/(2 - par)
        tau[par == 2] <- 1 - trigamma(2)
    } else if (family == 7 || family == 17) {
        theta <- par
        delta <- par2
        tau <- 1 - 2/(delta * (theta + 2))
    } else if (family == 8 || family == 18) {
        theta <- par
        delta <- par2
        kt <- function(t, th, de) {
            -log(-(1 - t)^th + 1) * (1 - t - (1 - t)^(-th) + (1 - t)^(-th) * t)/(de * th)
        }
        tau <- 1 + 4 * mapply(function(theta, delta) {
            integrate(function(t) {
                kt(t, th = theta, de = delta)
            }, 0, 1)$value
        }, theta, delta)
    } else if (family == 9 || family == 19) {
        theta <- par
        delta <- par2
        # tau=1-2/(delta*(2-theta))+4/(theta^2*delta)*gamma(delta+2)*gamma((2-2*theta)/(theta)+1)/gamma(delta+3+(2-2*theta)/(theta))
        kt <- function(t, th, de) {
            ((1 - (1 - t)^th)^-de - 1)/(-th * de * (1 - t)^(th - 1) * (1 - (1 - t)^th)^(-de - 1))
        }
        tau <- 1 + 4 * mapply(function(theta, delta) {
            integrate(function(t)  kt(t, th = theta, de = delta), 0, 1)$value
        }, theta, delta)
        # kt <- function(t) { 1/( (1-t)^(par-1) ) } kt2 <- function(t) { 1-t } kt3 <-
        # function(t) { 1/( (1-t)^(par-1)*(1-(1-t)^par)^(-par2-1) ) } tau <-
        # 1-4/par/par2*(integrate(kt,0,1)$value-integrate(kt2,0,1)$value-integrate(kt3,0,1)$value)
    } else if (family == 10 || family == 20) {
        theta <- par
        delta <- par2
        kt <- function(t, th, de) {
            -log(((1 - t * de)^th - 1)/((1 - de)^th - 1)) * (1 - t * de - (1 - t * de)^(-th) + (1 - t * de)^(-th) * t * de)/(th * de)
        }
        tau <- 1 + 4 * mapply(function(theta, delta) {
            integrate(function(t) {
                kt(t, th = theta, de = delta)
            }, 0, 1)$value
        }, theta, delta)
    } else if (family == 23 || family == 33) {
        tau <- par/(-par + 2)
    } else if (family == 24 || family == 34) {
        tau <- -1 - 1/par
    } else if (family == 26 || family == 36) {
        # tau <- -1-4/par^2*integrate(function(x) log(x)*x*(1-x)^(2*(1+par)/-par), 0,
        # 1)$value
        theta <- -par
        param1 <- 2/theta + 1
        tem <- digamma(2) - digamma(param1)
        tau <- 1 + tem * 2/(2 - theta)
        tau[theta == 2] <- 1 - trigamma(2)
        tau <- -tau
    } else if (family == 27 || family == 37) {
        theta <- -par
        delta <- -par2
        tau <- 1 - 2/(delta * (theta + 2))
        tau <- -tau
    } else if (family == 28 || family == 38) {
        theta <- -par
        delta <- -par2
        kt <- function(t, th, de) {
            -log(-(1 - t)^th + 1) * (1 - t - (1 - t)^(-th) + (1 - t)^(-th) * t)/(de * th)
        }
        tau <- 1 + 4 * mapply(function(theta, delta) {
            integrate(function(t) {
                kt(t, th = theta, de = delta)
            }, 0, 1)$value
        }, theta, delta)
        tau <- -tau
    } else if (family == 29 || family == 39) {
        theta <- -par
        delta <- -par2
        # tau <-
        # 1-2/(delta*(2-theta))+4/(theta^2*delta)*gamma(delta+2)*gamma((2-2*theta)/(theta)+1)/gamma(delta+3+(2-2*theta)/(theta))
        kt <- function(t, th, de) {
            ((1 - (1 - t)^th)^(-de) - 1)/(-th * de * (1 - t)^(th - 1) * (1 - (1 - t)^th)^(-de - 1))
        }
        tau <- 1 + 4 * mapply(function(theta, delta) {
            integrate(function(t) {
                kt(t, th = theta, de = delta)
            }, 0, 1)$value
        }, theta, delta)
        tau <- -tau
    } else if (family == 30 || family == 40) {
        theta <- -par
        delta <- -par2
        kt <- function(t, th, de) {
            -log(((1 - t * de)^th - 1)/((1 - de)^th - 1)) * (1 - t * de - (1 - t * de)^(-th) + (1 - t * de)^(-th) * t * de)/(th * de)
        }
        tau <- 1 + 4 * mapply(function(theta, delta) {
            integrate(function(t) {
                kt(t, th = theta, de = delta)
            }, 0, 1)$value
        }, theta, delta)
        tau <- -tau
    } else if (family == 41 || family == 51) {
        de <- par
        ln2 <- log(2)
        tem <- (2 - 2 * de) * ln2 + lgamma(2 * de) - 2 * lgamma(1 + de)
        tau <- 1 - de * exp(tem)
    } else if (family == 61 || family == 71) {
        de <- -par
        ln2 <- log(2)
        tem <- (2 - 2 * de) * ln2 + lgamma(2 * de) - 2 * lgamma(1 + de)
        tau <- 1 - de * exp(tem)
        tau <- -tau
    } else if (family == 42) {
        tau <- (75 * par2 - par2^2 + par * (25 - par2))/450
    } else if (family == 104 || family == 114 || family == 204 || family == 214) {
        par3 <- 1
        tau_int <- function(t, th, de) {
            Afunc <- .C("Tawn2",
                        as.double(t), 
                        as.integer(length(t)), 
                        as.double(th), 
                        as.double(de), 
                        as.double(1), 
                        as.double(rep(0, length(t))), 
                        PACKAGE = "VineCopula")[[6]]
            Afunc2Deriv <- .C("d2Tawn", 
                              as.double(t), 
                              as.integer(length(t)), 
                              as.double(th), 
                              as.double(de), 
                              as.double(1),
                              as.double(rep(0, length(t))),
                              PACKAGE = "VineCopula")[[6]]
            (t * (1 - t)) * Afunc2Deriv/Afunc
        }
        tau <- mapply(function(par, par2) {
            integrate(function(t) {
                tau_int(t, th = par, de = par2)
            }, 0, 1)$value
        }, par, par2)
    } else if (family == 124 || family == 134 || family == 224 || family == 234) {
        par3 <- 1
        tau_int <- function(t, th, de) {
            Afunc <- .C("Tawn2", 
                        as.double(t),
                        as.integer(length(t)),
                        as.double(-th), 
                        as.double(de),
                        as.double(1),
                        as.double(rep(0, length(t))),
                        PACKAGE = "VineCopula")[[6]]
            Afunc2Deriv <- .C("d2Tawn",
                              as.double(t),
                              as.integer(length(t)), 
                              as.double(-th), 
                              as.double(de), 
                              as.double(1),
                              as.double(rep(0, length(t))), 
                              PACKAGE = "VineCopula")[[6]]
            (t * (1 - t)) * Afunc2Deriv/Afunc
        }
        tau <- mapply(function(par, par2) {
            integrate(function(t) {
                tau_int(t, th = par, de = par2)
            }, 0, 1)$value
        }, par, par2)
        tau <- -tau
    }
    
    return(tau)
}
