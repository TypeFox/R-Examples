BiCopPar2TailDep <- function(family, par, par2 = 0) {
  if (!(family %in% c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 13, 14, 16, 17, 18, 19, 20, 23, 24, 26, 27, 28, 
    29, 30, 33, 34, 36, 37, 38, 39, 40))) 
    stop("Copula family not implemented.")
  if (c(2, 7, 8, 9, 10, 17, 18, 19, 20, 27, 28, 29, 30, 37, 38, 39, 40) %in% family && par2 == 0) 
    stop("For t-, BB1, BB6, BB7 and BB8 copulas, 'par2' must be set.")
  if (c(1, 3, 4, 5, 6, 13, 14, 16, 23, 24, 26, 33, 34, 36) %in% family && length(par) < 1) 
    stop("'par' not set.")
  
  if ((family == 1 || family == 2) && abs(par[1]) >= 1) 
    stop("The parameter of the Gaussian and t-copula has to be in the interval (-1,1).")
  if (family == 2 && par2 <= 2) 
    stop("The degrees of freedom parameter of the t-copula has to be larger than 2.")
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
  
  if (family == 0 | family == 1 | family == 5 | family %in% c(23, 24, 26, 27, 28, 29, 30, 33, 34, 36, 37, 
    38, 39, 40)) {
    lower <- 0
    upper <- 0
  } else if (family == 2) {
    lower <- 2 * pt((-sqrt(par2 + 1) * sqrt((1 - par)/(1 + par))), df = par2 + 1)
    upper <- lower
  } else if (family == 3) {
    lower <- 2^(-1/par)
    upper <- 0
  } else if (family == 4 | family == 6) {
    lower <- 0
    upper <- 2 - 2^(1/par)
  } else if (family == 7) {
    lower <- 2^(-1/(par * par2))
    upper <- 2 - 2^(1/par2)
  } else if (family == 8) {
    lower <- 0
    upper <- 2 - 2^(1/(par * par2))
  } else if (family == 9) {
    lower <- 2^(-1/par2)
    upper <- 2 - 2^(1/par)
  } else if (family == 10) {
    lower <- 0
    if (par2 == 1) 
      upper <- 2 - 2^(1/par) else upper <- 0
  } else if (family == 13) {
    lower <- 0
    upper <- 2^(-1/par)
  } else if (family == 14 | family == 16) {
    lower <- 2 - 2^(1/par)
    upper <- 0
  } else if (family == 17) {
    lower <- 2 - 2^(1/par2)
    upper <- 2^(-1/par2)
  } else if (family == 18) {
    lower <- 2 - 2^(1/(par * par2))
    upper <- 0
  } else if (family == 19) {
    lower <- 2 - 2^(1/par)
    upper <- 2^(-1/par2)
  } else if (family == 20) {
    if (par2 == 1) 
      lower <- 2 - 2^(1/par) else lower <- 0
    
    upper <- 0
  }
  
  return(list(lower = lower, upper = upper))
} 
