CDVineClarkeTest <- function(data, Model1.order = 1:dim(data)[2], Model2.order = 1:dim(data)[2],
                             Model1.family, Model2.family, Model1.par, Model2.par,
                             Model1.par2 = rep(0, dim(data)[2] * (dim(data)[2] - 1)/2),
                             Model2.par2 = rep(0, dim(data)[2] * (dim(data)[2] - 1)/2),
                             Model1.type, Model2.type) {
  data <- as.matrix(data)
  N <- dim(data)[1]
  d <- dim(data)[2]
  if (any(data > 1) || any(data < 0)) 
    stop("Data has be in the interval [0,1].")
  
  if (Model1.type == "CVine") 
    Model1.type <- 1 else if (Model1.type == "DVine") 
    Model1.type <- 2
  if (Model1.type != 1 & Model1.type != 2) 
    stop("Model1.type: Vine model not implemented.")
  
  if (Model2.type == "CVine") 
    Model2.type <- 1 else if (Model2.type == "DVine") 
    Model2.type <- 2
  if (Model2.type != 1 & Model2.type != 2) 
    stop("Model2.type: Vine model not implemented.")
  
  # Sicherheitsabfragen
  
  if (d < 2) 
    stop("Dimension has to be at least 2.")
  if (N < 2) 
    stop("Number of observations has to be at least 2.")
  
  
  if (length(Model1.family) != d * (d - 1)/2) 
    stop("Number of copula families incorrect in model 1.")
  if (length(Model2.family) != d * (d - 1)/2) 
    stop("Number of copula families incorrect in model 2.")
  if (length(Model1.par) != d * (d - 1)/2) 
    stop("Number of copula parameters incorrect in model 1.")
  if (length(Model2.par) != d * (d - 1)/2) 
    stop("Number of copula parameters incorrect in model 2.")
  if (length(Model1.par2) != d * (d - 1)/2) 
    stop("Number of second copula parameters incorrect in model 1.")
  if (length(Model2.par2) != d * (d - 1)/2) 
    stop("Number of second copula parameters incorrect in model 2.")
  
  if (length(Model1.order) != d) 
    stop("Length of order vector incorrect in model 1.")
  if (length(Model2.order) != d) 
    stop("Length of order vector incorrect in model 2.")
  
  if (max(Model1.order) > d || min(Model1.order) < 1) 
    stop("Order vector of model 1 is incorrect.")
  if (max(Model2.order) > d || min(Model2.order) < 1) 
    stop("Order vector of model 2 is incorrect.")
  
  for (i in 1:(d * (d - 1)/2)) {
    if (!(Model1.family[i] %in% c(0, 1:10, 13, 14, 16:20, 23, 24, 26:30, 33, 34, 36:40))) 
      stop("Copula family not implemented (in Model1.family).")
    if ((Model1.family[i] == 1 || Model1.family[i] == 2) && abs(Model1.par[i]) >= 1) 
      stop("The parameter of the Gaussian and t-copula has to be in the interval (-1,1).")
    if (Model1.family[i] == 2 && Model1.par2[i] <= 2) 
      stop("The degrees of freedom parameter of the t-copula has to be larger than 2.")
    if ((Model1.family[i] == 3 || Model1.family[i] == 13) && Model1.par[i] <= 0) 
      stop("The parameter of the Clayton copula has to be positive.")
    if ((Model1.family[i] == 4 || Model1.family[i] == 14) && Model1.par[i] < 1) 
      stop("The parameter of the Gumbel copula has to be in the interval [1,oo).")
    if ((Model1.family[i] == 6 || Model1.family[i] == 16) && Model1.par[i] <= 1) 
      stop("The parameter of the Joe copula has to be in the interval (1,oo).")
    if (Model1.family[i] == 5 && Model1.par[i] == 0) 
      stop("The parameter of the Frank copula has to be unequal to 0.")
    if ((Model1.family[i] == 7 || Model1.family[i] == 17) && Model1.par[i] <= 0) 
      stop("The first parameter of the BB1 copula has to be positive.")
    if ((Model1.family[i] == 7 || Model1.family[i] == 17) && Model1.par2[i] < 1) 
      stop("The second parameter of the BB1 copula has to be in the interval [1,oo).")
    if ((Model1.family[i] == 8 || Model1.family[i] == 18) && Model1.par[i] <= 0) 
      stop("The first parameter of the BB6 copula has to be in the interval [1,oo).")
    if ((Model1.family[i] == 8 || Model1.family[i] == 18) && Model1.par2[i] < 1) 
      stop("The second parameter of the BB6 copula has to be in the interval [1,oo).")
    if ((Model1.family[i] == 9 || Model1.family[i] == 19) && Model1.par[i] < 1) 
      stop("The first parameter of the BB7 copula has to be in the interval [1,oo).")
    if ((Model1.family[i] == 9 || Model1.family[i] == 19) && Model1.par2[i] <= 0) 
      stop("The second parameter of the BB7 copula has to be positive.")
    if ((Model1.family[i] == 10 || Model1.family[i] == 20) && Model1.par[i] < 1) 
      stop("The first parameter of the BB8 copula has to be in the interval [1,oo).")
    if ((Model1.family[i] == 10 || Model1.family[i] == 20) && (Model1.par2[i] <= 0 || Model1.par2[i] > 
      1)) 
      stop("The second parameter of the BB8 copula has to be in the interval (0,1].")
    if ((Model1.family[i] == 23 || Model1.family[i] == 33) && Model1.par[i] >= 0) 
      stop("The parameter of the rotated Clayton copula has to be negative.")
    if ((Model1.family[i] == 24 || Model1.family[i] == 34) && Model1.par[i] > -1) 
      stop("The parameter of the rotated Gumbel copula has to be in the interval (-oo,-1].")
    if ((Model1.family[i] == 26 || Model1.family[i] == 36) && Model1.par[i] >= -1) 
      stop("The parameter of the rotated Joe copula has to be in the interval (-oo,-1).")
    if ((Model1.family[i] == 27 || Model1.family[i] == 37) && Model1.par[i] >= 0) 
      stop("The first parameter of the rotated BB1 copula has to be negative.")
    if ((Model1.family[i] == 27 || Model1.family[i] == 37) && Model1.par2[i] > -1) 
      stop("The second parameter of the rotated BB1 copula has to be in the interval (-oo,-1].")
    if ((Model1.family[i] == 28 || Model1.family[i] == 38) && Model1.par[i] >= 0) 
      stop("The first parameter of the rotated BB6 copula has to be in the interval (-oo,-1].")
    if ((Model1.family[i] == 28 || Model1.family[i] == 38) && Model1.par2[i] > -1) 
      stop("The second parameter of the rotated BB6 copula has to be in the interval (-oo,-1].")
    if ((Model1.family[i] == 29 || Model1.family[i] == 39) && Model1.par[i] > -1) 
      stop("The first parameter of the rotated BB7 copula has to be in the interval (-oo,-1].")
    if ((Model1.family[i] == 29 || Model1.family[i] == 39) && Model1.par2[i] >= 0) 
      stop("The second parameter of the rotated BB7 copula has to be negative.")
    if ((Model1.family[i] == 30 || Model1.family[i] == 40) && Model1.par[i] > -1) 
      stop("The first parameter of the rotated BB8 copula has to be in the interval (-oo,-1].")
    if ((Model1.family[i] == 30 || Model1.family[i] == 40) && (Model1.par2[i] >= 0 || Model1.par2[i] < 
      (-1))) 
      stop("The second parameter of the rotated BB8 copula has to be in the interval [-1,0).")
  }
  for (i in 1:(d * (d - 1)/2)) {
    if (!(Model2.family[i] %in% c(0, 1:10, 13, 14, 16:20, 23, 24, 26:30, 33, 34, 36:40))) 
      stop("Copula family not implemented (in Model1.family).")
    if ((Model2.family[i] == 1 || Model2.family[i] == 2) && abs(Model2.par[i]) >= 1) 
      stop("The parameter of the Gaussian and t-copula has to be in the interval (-1,1).")
    if (Model2.family[i] == 2 && Model2.par2[i] <= 2) 
      stop("The degrees of freedom parameter of the t-copula has to be larger than 2.")
    if ((Model2.family[i] == 3 || Model2.family[i] == 13) && Model2.par[i] <= 0) 
      stop("The parameter of the Clayton copula has to be positive.")
    if ((Model2.family[i] == 4 || Model2.family[i] == 14) && Model2.par[i] < 1) 
      stop("The parameter of the Gumbel copula has to be in the interval [1,oo).")
    if ((Model2.family[i] == 6 || Model2.family[i] == 16) && Model2.par[i] <= 1) 
      stop("The parameter of the Joe copula has to be in the interval (1,oo).")
    if (Model2.family[i] == 5 && Model2.par[i] == 0) 
      stop("The parameter of the Frank copula has to be unequal to 0.")
    if ((Model2.family[i] == 7 || Model2.family[i] == 17) && Model2.par[i] <= 0) 
      stop("The first parameter of the BB1 copula has to be positive.")
    if ((Model2.family[i] == 7 || Model2.family[i] == 17) && Model2.par2[i] < 1) 
      stop("The second parameter of the BB1 copula has to be in the interval [1,oo).")
    if ((Model2.family[i] == 8 || Model2.family[i] == 18) && Model2.par[i] <= 0) 
      stop("The first parameter of the BB6 copula has to be in the interval [1,oo).")
    if ((Model2.family[i] == 8 || Model2.family[i] == 18) && Model2.par2[i] < 1) 
      stop("The second parameter of the BB6 copula has to be in the interval [1,oo).")
    if ((Model2.family[i] == 9 || Model2.family[i] == 19) && Model2.par[i] < 1) 
      stop("The first parameter of the BB7 copula has to be in the interval [1,oo).")
    if ((Model2.family[i] == 9 || Model2.family[i] == 19) && Model2.par2[i] <= 0) 
      stop("The second parameter of the BB7 copula has to be positive.")
    if ((Model2.family[i] == 10 || Model2.family[i] == 20) && Model2.par[i] < 1) 
      stop("The first parameter of the BB8 copula has to be in the interval [1,oo).")
    if ((Model2.family[i] == 10 || Model2.family[i] == 20) && (Model2.par2[i] <= 0 || Model2.par2[i] > 
      1)) 
      stop("The second parameter of the BB8 copula has to be in the interval (0,1].")
    if ((Model2.family[i] == 23 || Model2.family[i] == 33) && Model2.par[i] >= 0) 
      stop("The parameter of the rotated Clayton copula has to be negative.")
    if ((Model2.family[i] == 24 || Model2.family[i] == 34) && Model2.par[i] > -1) 
      stop("The parameter of the rotated Gumbel copula has to be in the interval (-oo,-1].")
    if ((Model2.family[i] == 26 || Model2.family[i] == 36) && Model2.par[i] >= -1) 
      stop("The parameter of the rotated Joe copula has to be in the interval (-oo,-1).")
    if ((Model2.family[i] == 27 || Model2.family[i] == 37) && Model2.par[i] >= 0) 
      stop("The first parameter of the rotated BB1 copula has to be negative.")
    if ((Model2.family[i] == 27 || Model2.family[i] == 37) && Model2.par2[i] > -1) 
      stop("The second parameter of the rotated BB1 copula has to be in the interval (-oo,-1].")
    if ((Model2.family[i] == 28 || Model2.family[i] == 38) && Model2.par[i] >= 0) 
      stop("The first parameter of the rotated BB6 copula has to be in the interval (-oo,-1].")
    if ((Model2.family[i] == 28 || Model2.family[i] == 38) && Model2.par2[i] > -1) 
      stop("The second parameter of the rotated BB6 copula has to be in the interval (-oo,-1].")
    if ((Model2.family[i] == 29 || Model2.family[i] == 39) && Model2.par[i] > -1) 
      stop("The first parameter of the rotated BB7 copula has to be in the interval (-oo,-1].")
    if ((Model2.family[i] == 29 || Model2.family[i] == 39) && Model2.par2[i] >= 0) 
      stop("The second parameter of the rotated BB7 copula has to be negative.")
    if ((Model2.family[i] == 30 || Model2.family[i] == 40) && Model2.par[i] > -1) 
      stop("The first parameter of the rotated BB8 copula has to be in the interval (-oo,-1].")
    if ((Model2.family[i] == 30 || Model2.family[i] == 40) && (Model2.par2[i] >= 0 || Model2.par2[i] < 
      (-1))) 
      stop("The second parameter of the rotated BB8 copula has to be in the interval [-1,0).")
  }
  
  
  
  Model1.ll <- numeric(N)
  Model2.ll <- numeric(N)
  
  for (i in 1:N) {
    Model1.ll[i] <- LogLikelihood(Model1.family, data[i, Model1.order], Model1.par, Model1.par2, Model1.type)
    Model2.ll[i] <- LogLikelihood(Model2.family, data[i, Model2.order], Model2.par, Model2.par2, Model2.type)
  }
  
  # model1.ll<-sapply(data,function (x) LogLikelihood(Model1.family, x, Model1.par, Model1.par2,
  # Model1.type))
  
  
  anz.1 <- sum(Model1.family >= 1, na.rm = TRUE) + sum(Model1.family %in% c(2, 7:10, 17:20, 27:30, 37:40), 
    na.rm = TRUE)
  anz.2 <- sum(Model2.family >= 1, na.rm = TRUE) + sum(Model2.family %in% c(2, 7:10, 17:20, 27:30, 37:40), 
    na.rm = TRUE)
  
  B <- sum(Model1.ll - Model2.ll > 0)
  B.Schwarz <- sum(Model1.ll - Model2.ll - (anz.1 - anz.2) * log(N)/(2 * N) > 0)
  B.Akaike <- sum(Model1.ll - Model2.ll - (anz.1 - anz.2)/N > 0)
  
  if (B == 0 | B == N/2) 
    p <- 1 else p <- 2 * min(pbinom(B, N, 0.5), 1 - pbinom(B - 1, N, 0.5))
  
  if (B.Schwarz == 0 | B.Schwarz == N/2) 
    p.Schwarz <- 1 else p.Schwarz <- 2 * min(pbinom(B.Schwarz, N, 0.5), 1 - pbinom(B.Schwarz - 1, N, 0.5))
  
  if (B.Akaike == 0 | B.Akaike == N/2) 
    p.Akaike <- 1 else p.Akaike <- 2 * min(pbinom(B.Akaike, N, 0.5), 1 - pbinom(B.Akaike - 1, N, 0.5))
  
  Vuong <- list(statistic = B, statistic.Akaike = B.Akaike, statistic.Schwarz = B.Schwarz, p.value = p, 
    p.value.Akaike = p.Akaike, p.value.Schwarz = p.Schwarz)
  
  return(Vuong)
  
}


#################################################
# LogLikelihood  				#
# Input:					#
# data		data vector			#
# cop		copula family vector		#
# par		first parameter vector		#
# nu		second parameter vector		#
# type		vine type			#
# Output:					#
# loglik	log-likelihood			#
#################################################

LogLikelihood <- function(cop, data, par, nu = rep(0, length(data)), type) {
  loglik <- 0
  d <- length(data)
  T <- 1
  ll <- rep(0, d * (d - 1)/2)
  if (type == 1) 
    vv <- rep(0, T * (((d - 1) * d)/2 - 1)) else vv <- rep(0, T * (d - 1) * (d - 2))
  out <- .C("VineLogLikm", as.integer(1), as.integer(d), as.integer(type), as.integer(cop), as.double(c(par, 
    nu)), as.double(data), as.double(0), as.double(ll), as.double(vv), PACKAGE = "CDVine")
  loglik <- -out[[7]]
  if (loglik %in% c(NA, NaN, -Inf)) 
    loglik <- -1e+10
  
  return(loglik)
} 
