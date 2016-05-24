CDVineLogLik <- function(data, family, par, par2 = rep(0, dim(data)[2] * (dim(data)[2] - 1)/2), type) {
  if (is.vector(data)) {
    data <- t(as.matrix(data))
  } else {
    data <- as.matrix(data)
  }
  if (any(data > 1) || any(data < 0)) 
    stop("Data has be in the interval [0,1].")
  if (type == "CVine") 
    type <- 1 else if (type == "DVine") 
    type <- 2
  if (type != 1 & type != 2) 
    stop("Vine model not implemented.")
  
  d <- dim(data)[2]
  T <- dim(data)[1]
  
  # Sicherheitsabfragen
  if (d < 2) 
    stop("Dimension has to be at least 2.")
  
  if (length(family) != d * (d - 1)/2) 
    stop("Number of copula families incorrect.")
  if (length(par) != d * (d - 1)/2) 
    stop("Number of copula parameters incorrect.")
  if (length(par2) != d * (d - 1)/2) 
    stop("Number of second copula parameters incorrect.")
  
  for (i in 1:(d * (d - 1)/2)) {
    if (!(family[i] %in% c(0, 1:10, 13, 14, 16:20, 23, 24, 26:30, 33, 34, 36:40))) 
      stop("Copula family not implemented.")
    # Parameterbereiche abfragen
    if ((family[i] == 1 || family[i] == 2) && abs(par[i]) >= 1) 
      stop("The parameter of the Gaussian and t-copula has to be in the interval (-1,1).")
    if (family[i] == 2 && par2[i] <= 2) 
      stop("The degrees of freedom parameter of the t-copula has to be larger than 2.")
    if ((family[i] == 3 || family[i] == 13) && par[i] <= 0) 
      stop("The parameter of the Clayton copula has to be positive.")
    if ((family[i] == 4 || family[i] == 14) && par[i] < 1) 
      stop("The parameter of the Gumbel copula has to be in the interval [1,oo).")
    if ((family[i] == 6 || family[i] == 16) && par[i] <= 1) 
      stop("The copula parameter of the Joe copula has to be in the interval (1,oo).")
    if (family[i] == 5 && par[i] == 0) 
      stop("The parameter of the Frank copula has to be unequal to 0.")
    if ((family[i] == 7 || family[i] == 17) && par[i] <= 0) 
      stop("The first parameter of the BB1 copula has to be positive.")
    if ((family[i] == 7 || family[i] == 17) && par2[i] < 1) 
      stop("The second parameter of the BB1 copula has to be in the interval [1,oo).")
    if ((family[i] == 8 || family[i] == 18) && par[i] <= 0) 
      stop("The first parameter of the BB6 copula has to be in the interval [1,oo).")
    if ((family[i] == 8 || family[i] == 18) && par2[i] < 1) 
      stop("The second parameter of the BB6 copula has to be in the interval [1,oo).")
    if ((family[i] == 9 || family[i] == 19) && par[i] < 1) 
      stop("The first parameter of the BB7 copula has to be in the interval [1,oo).")
    if ((family[i] == 9 || family[i] == 19) && par2[i] <= 0) 
      stop("The second parameter of the BB7 copula has to be positive.")
    if ((family[i] == 10 || family[i] == 20) && par[i] < 1) 
      stop("The first parameter of the BB8 copula has to be in the interval [1,oo).")
    if ((family[i] == 10 || family[i] == 20) && (par2[i] <= 0 || par2[i] > 1)) 
      stop("The second parameter of the BB8 copula has to be in the interval (0,1].")
    if ((family[i] == 23 || family[i] == 33) && par[i] >= 0) 
      stop("The parameter of the rotated Clayton copula has to be negative.")
    if ((family[i] == 24 || family[i] == 34) && par[i] > -1) 
      stop("The parameter of the rotated Gumbel copula has to be in the interval (-oo,-1].")
    if ((family[i] == 26 || family[i] == 36) && par[i] >= -1) 
      stop("The parameter of the rotated Joe copula has to be in the interval (-oo,-1).")
    if ((family[i] == 27 || family[i] == 37) && par[i] >= 0) 
      stop("The first parameter of the rotated BB1 copula has to be negative.")
    if ((family[i] == 27 || family[i] == 37) && par2[i] > -1) 
      stop("The second parameter of the rotated BB1 copula has to be in the interval (-oo,-1].")
    if ((family[i] == 28 || family[i] == 38) && par[i] >= 0) 
      stop("The first parameter of the rotated BB6 copula has to be in the interval (-oo,-1].")
    if ((family[i] == 28 || family[i] == 38) && par2[i] > -1) 
      stop("The second parameter of the rotated BB6 copula has to be in the interval (-oo,-1].")
    if ((family[i] == 29 || family[i] == 39) && par[i] > -1) 
      stop("The first parameter of the rotated BB7 copula has to be in the interval (-oo,-1].")
    if ((family[i] == 29 || family[i] == 39) && par2[i] >= 0) 
      stop("The second parameter of the rotated BB7 copula has to be negative.")
    if ((family[i] == 30 || family[i] == 40) && par[i] > -1) 
      stop("The first parameter of the rotated BB8 copula has to be in the interval (-oo,-1].")
    if ((family[i] == 30 || family[i] == 40) && (par2[i] >= 0 || par2[i] < (-1))) 
      stop("The second parameter of the rotated BB8 copula has to be in the interval [-1,0).")
  }
  
  ll <- rep(0, d * (d - 1)/2)
  if (d > 2) {
    if (type == 1) 
      vv <- rep(0, T * (((d - 1) * d)/2 - 1)) else vv <- rep(0, T * (d - 1) * (d - 2))
  } else {
    vv <- 0
  }
  w1 <- family
  
  th <- c(par, par2)
  
  if (d > 2) {
    
    out <- .C("VineLogLikm", as.integer(T), as.integer(d), as.integer(type), as.integer(w1), as.double(th), 
      as.double(data), as.double(0), as.double(ll), as.double(vv), PACKAGE = "CDVine")
    ll <- out[[8]]
    ll[ll %in% c(NA, NaN, -Inf)] <- 1e+10
    loglik <- -out[[7]]
    
  } else {
    
    ll <- log(BiCopPDF(data[, 1], data[, 2], family, par, par2))
    loglik <- sum(ll)
    
  }
  
  
  if (loglik %in% c(NA, NaN, -Inf)) 
    loglik <- -1e+10
  if (d > 2) {
    vv <- out[[9]]
    vv[vv > 0.99999999] <- 0.99999999
    vv[vv < 1e-08] <- 1e-08
  }
  res <- list(loglik = loglik, ll = ll, vv = vv)
  return(res)
} 
