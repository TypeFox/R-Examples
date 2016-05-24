BiCopCDF <- function(u1, u2, family, par, par2 = 0) {
  if (is.null(u1) == TRUE || is.null(u2) == TRUE) 
    stop("u1 and/or u2 are not set or have length zero.")
  if (any(u1 > 1) || any(u1 < 0)) 
    stop("Data has be in the interval [0,1].")
  if (any(u2 > 1) || any(u2 < 0)) 
    stop("Data has be in the interval [0,1].")
  if (length(u1) != length(u2)) 
    stop("Lengths of 'u1' and 'u2' do not match.")
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
  
  res <- rep(NA, length(u1))
  
  if (family == 0) {
    res <- u1 * u2
  } else if (family == 1) {
    cdf <- function(u, v) pmvnorm(upper = c(qnorm(u), qnorm(v)), corr = matrix(c(1, par, par, 1), 2, 
      2))
    res <- mapply(cdf, u1, u2, SIMPLIFY = TRUE)
  } else if (family == 2) {
    cdf <- function(u, v) pmvt(upper = c(qt(u, df = par2), qt(v, df = par2)), corr = matrix(c(1, par, 
      par, 1), 2, 2), df = par2)
    res <- mapply(cdf, u1, u2, SIMPLIFY = TRUE)
  } else if (family %in% c(3:10)) {
    res <- .C("archCDF", as.double(u1), as.double(u2), as.integer(length(u1)), as.double(c(par, par2)), 
      as.integer(family), as.double(rep(0, length(u1))), PACKAGE = "CDVine")[[6]]
  } else if (family %in% c(13, 14, 16:20)) {
    res <- u1 + u2 - 1 + .C("archCDF", as.double(1 - u1), as.double(1 - u2), as.integer(length(u1)), 
      as.double(c(par, par2)), as.integer(family - 10), as.double(rep(0, length(u1))), PACKAGE = "CDVine")[[6]]
  } else if (family %in% c(23, 24, 26:30)) {
    res <- u2 - .C("archCDF", as.double(1 - u1), as.double(u2), as.integer(length(u1)), as.double(c(-par, 
      -par2)), as.integer(family - 20), as.double(rep(0, length(u1))), PACKAGE = "CDVine")[[6]]
  } else if (family %in% c(33, 34, 36:40)) {
    res <- u1 - .C("archCDF", as.double(u1), as.double(1 - u2), as.integer(length(u1)), as.double(c(-par, 
      -par2)), as.integer(family - 30), as.double(rep(0, length(u1))), PACKAGE = "CDVine")[[6]]
  }
  
  return(res)
} 
