BiCopPar2Tau <- function(family, par, par2 = 0) {
  if (!(family %in% c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 13, 14, 16, 17, 18, 19, 20, 23, 24, 26, 27, 28, 
    29, 30, 33, 34, 36, 37, 38, 39, 40))) 
    stop("Copula family not implemented.")
  if (family %in% c(2, 7, 8, 9, 10, 17, 18, 19, 20, 27, 28, 29, 30, 37, 38, 39, 40) && par2 == 0) 
    stop("For t-, BB1, BB6, BB7 and BB8 copulas, 'par2' must be set.")
  if (family %in% c(1, 3, 4, 5, 6, 13, 14, 16, 23, 24, 26, 33, 34, 36) && length(par) < 1) 
    stop("'par' not set.")
  
  if ((family == 1 || family == 2) && abs(par[1]) >= 1) 
    stop("The parameter of the Gaussian and t-copula has to be in the interval (-1,1).")
  # if(family==2 && par2<=1) stop('The degrees of freedom parameter of the t-copula has to be larger than
  # 1.')
  if ((family == 3 || family == 13) && par <= 0) 
    stop("The parameter of the Clayton copula has to be positive.")
  if ((family == 4 || family == 14) && par < 1) 
    stop("The parameter of the Gumbel copula has to be in the interval [1,oo).")
  if ((family == 6 || family == 16) && par <= 1) 
    stop("The parameter of the Joe copula has to be in the interval (1,oo).")
  if (family == 5 && par == 0) 
    stop("The parameter of the Frank copula has to be unequal to 0.")
  if ((family == 7 || family == 17) && par <= 0) 
    stop("The first parameter of the BB1 copula has to be positive.")
  if ((family == 7 || family == 17) && par2 < 1) 
    stop("The second parameter of the BB1 copula has to be in the interval [1,oo).")
  if ((family == 8 || family == 18) && par <= 0) 
    stop("The first parameter of the BB6 copula has to be in the interval [1,oo).")
  if ((family == 8 || family == 18) && par2 < 1) 
    stop("The second parameter of the BB6 copula has to be in the interval [1,oo).")
  if ((family == 9 || family == 19) && par < 1) 
    stop("The first parameter of the BB7 copula has to be in the interval [1,oo).")
  if ((family == 9 || family == 19) && par2 <= 0) 
    stop("The second parameter of the BB7 copula has to be positive.")
  if ((family == 10 || family == 20) && par < 1) 
    stop("The first parameter of the BB8 copula has to be in the interval [1,oo).")
  if ((family == 10 || family == 20) && (par2 <= 0 || par2 > 1)) 
    stop("The second parameter of the BB8 copula has to be in the interval (0,1].")
  if ((family == 23 || family == 33) && par >= 0) 
    stop("The parameter of the rotated Clayton copula has to be negative.")
  if ((family == 24 || family == 34) && par > -1) 
    stop("The parameter of the rotated Gumbel copula has to be in the interval (-oo,-1].")
  if ((family == 26 || family == 36) && par >= -1) 
    stop("The parameter of the rotated Joe copula has to be in the interval (-oo,-1).")
  if ((family == 27 || family == 37) && par >= 0) 
    stop("The first parameter of the rotated BB1 copula has to be negative.")
  if ((family == 27 || family == 37) && par2 > -1) 
    stop("The second parameter of the rotated BB1 copula has to be in the interval (-oo,-1].")
  if ((family == 28 || family == 38) && par >= 0) 
    stop("The first parameter of the rotated BB6 copula has to be in the interval (-oo,-1].")
  if ((family == 28 || family == 38) && par2 > -1) 
    stop("The second parameter of the rotated BB6 copula has to be in the interval (-oo,-1].")
  if ((family == 29 || family == 39) && par > -1) 
    stop("The first parameter of the rotated BB7 copula has to be in the interval (-oo,-1].")
  if ((family == 29 || family == 39) && par2 >= 0) 
    stop("The second parameter of the rotated BB7 copula has to be negative.")
  if ((family == 30 || family == 40) && par > -1) 
    stop("The first parameter of the rotated BB8 copula has to be in the interval (-oo,-1].")
  if ((family == 30 || family == 40) && (par2 >= 0 || par2 < (-1))) 
    stop("The second parameter of the rotated BB8 copula has to be in the interval [-1,0).")
  
  if (family == 0) {
    tau <- 0
  } else if (family == 1 | family == 2) {
    tau <- 2/pi * asin(par)
  } else if (family == 3 || family == 13) {
    tau <- par/(par + 2)
  } else if (family == 4 || family == 14) {
    tau <- 1 - 1/par
  } else if (family == 5) {
    f <- function(x) {
      x/(exp(x) - 1)
    }
    if (par > 0) {
      tau <- 1 - 4/par + 4/par^2 * integrate(f, lower = 0, upper = par)$value
    } else {
      tau <- 1 - 4/par - 4/par^2 * integrate(f, lower = par, upper = 0)$value
    }
  } else if (family == 6 || family == 16) {
    # tau=1+4/par^2*integrate(function(x) log(x)*x*(1-x)^(2*(1-par)/par), 0, 1)$value
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
    kt <- function(t) {
      -log(-(1 - t)^theta + 1) * (1 - t - (1 - t)^(-theta) + (1 - t)^(-theta) * t)/(delta * theta)
    }
    tau <- 1 + 4 * integrate(kt, 0, 1)$value
  } else if (family == 9 || family == 19) {
    theta <- par
    delta <- par2
    # tau=1-2/(delta*(2-theta))+4/(theta^2*delta)*gamma(delta+2)*gamma((2-2*theta)/(theta)+1)/gamma(delta+3+(2-2*theta)/(theta))
    kt <- function(t) {
      ((1 - (1 - t)^par)^-par2 - 1)/(-par * par2 * (1 - t)^(par - 1) * (1 - (1 - t)^par)^(-par2 - 1))
    }
    tau <- 1 + 4 * integrate(kt, 0, 1)$value
  } else if (family == 10 || family == 20) {
    theta <- par
    delta <- par2
    kt <- function(t) {
      -log(((1 - t * delta)^theta - 1)/((1 - delta)^theta - 1)) * (1 - t * delta - (1 - t * delta)^(-theta) + 
        (1 - t * delta)^(-theta) * t * delta)/(theta * delta)
    }
    tau <- 1 + 4 * integrate(kt, 0, 1)$value
  } else if (family == 23 || family == 33) {
    tau <- par/(-par + 2)
  } else if (family == 24 || family == 34) {
    tau <- -1 - 1/par
  } else if (family == 26 || family == 36) {
    # tau=-1-4/par^2*integrate(function(x) log(x)*x*(1-x)^(2*(1+par)/-par), 0, 1)$value
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
    kt <- function(t) {
      -log(-(1 - t)^theta + 1) * (1 - t - (1 - t)^(-theta) + (1 - t)^(-theta) * t)/(delta * theta)
    }
    tau <- 1 + 4 * integrate(kt, 0, 1)$value
    tau <- -tau
  } else if (family == 29 || family == 39) {
    theta <- -par
    delta <- -par2
    # tau=1-2/(delta*(2-theta))+4/(theta^2*delta)*gamma(delta+2)*gamma((2-2*theta)/(theta)+1)/gamma(delta+3+(2-2*theta)/(theta))
    kt <- function(t) {
      ((1 - (1 - t)^par)^(-par2) - 1)/(-par * par2 * (1 - t)^(par - 1) * (1 - (1 - t)^par)^(par2 - 
        1))
    }
    tau <- 1 + 4 * integrate(kt, 0, 1)$value
    tau <- -tau
  } else if (family == 30 || family == 40) {
    theta <- -par
    delta <- -par2
    kt <- function(t) {
      -log(((1 - t * delta)^theta - 1)/((1 - delta)^theta - 1)) * (1 - t * delta - (1 - t * delta)^(-theta) + 
        (1 - t * delta)^(-theta) * t * delta)/(theta * delta)
    }
    tau <- 1 + 4 * integrate(kt, 0, 1)$value
    tau <- -tau
  }
  
  return(tau)
} 
