#################################################################
# gen  							#
# Input:							#
# u		data vector					#
# param		copula parameters				#
# copula	copula family (7=BB1, 8=BB6, 9=BB7, 10=BB8)	#
# Output:							#
# out		generator					#
#################################################################

gen <- function(u, param, copula) {
  # param == c(theta, delta)
  out <- numeric(length(u))
  
  if (copula == 7) {
    out <- (u^(-param[1]) - 1)^param[2]
  } else if (copula == 8) {
    out <- (-log(-(1 - u)^param[1] + 1))^param[2]
  } else if (copula == 9) {
    out <- (1 - (1 - u)^param[1])^(-param[2]) - 1
  } else if (copula == 10) {
    out <- -log((1 - (1 - param[2] * u)^param[1])/(1 - (1 - param[2])^param[1]))
  }
  
  return(out)
}


#################################################################
# gen.inv  						#
# Input:							#
# u		data vector					#
# param		copula parameters				#
# copula	copula family (7=BB1, 8=BB6, 9=BB7, 10=BB8)	#
# Output:							#
# out		inverse generator				#
#################################################################

gen.inv <- function(u, param, copula) {
  out <- numeric(length(u))
  
  if (copula == 7) {
    out <- (1 + u^(1/param[2]))^(-1/param[1])
  } else if (copula == 8) {
    out <- 1 - (1 - exp(-u^(1/param[2])))^(1/param[1])
  } else if (copula == 9) {
    out <- 1 - (1 - (1 + u)^(-1/param[2]))^(1/param[1])
  } else if (copula == 10) {
    out <- 1/param[2] * (1 - (1 - (1 - (1 - param[2])^param[1]) * exp(-u))^(1/param[1]))
  }
  
  return(out)
}


#################################################################
# gen.drv  						#
# Input:							#
# u		data vector					#
# param		copula parameters				#
# copula	copula family (7=BB1, 8=BB6, 9=BB7, 10=BB8)	#
# Output:							#
# out		First derivative of the generator		#
#################################################################

gen.drv <- function(u, param, copula) {
  out <- numeric(length(u))
  
  if (copula == 7) {
    out <- -prod(param) * (u^-(param[1]) - 1)^(param[2] - 1) * u^(-1 - param[1])
  } else if (copula == 8) {
    out <- (-log(-(1 - u)^param[1] + 1))^(param[2] - 1) * param[2] * (1 - u)^(param[1] - 1) * param[1]/((1 - 
      u)^param[1] - 1)
  } else if (copula == 9) {
    out <- -prod(param) * (1 - u)^(param[1] - 1) * (1 - (1 - u)^param[1])^(-1 - param[2])
  } else if (copula == 10) {
    out <- -prod(param) * ((1 - param[2] * u)^(param[1] - 1))/(1 - (1 - param[2] * u)^param[1])
  }
  
  return(out)
}


#################################################################
# gen.drv2  						#
# Input:							#
# u		data vector					#
# param		copula parameters				#
# copula	copula family (7=BB1, 8=BB6, 9=BB7, 10=BB8)	#
# Output:							#
# out		Second derivative of the generator		#
#################################################################

gen.drv2 <- function(u, param, copula) {
  out <- numeric(length(u))
  
  if (copula == 7) {
    out <- prod(param) * u^(-2 - param[1]) * (u^(-param[1]) - 1)^(param[2] - 2) * ((1 + prod(param)) * 
      u^(-param[1]) - param[1] - 1)
  } else if (copula == 8) {
    out <- (prod(param) * ((-log(-(1 - u)^param[1] + 1))^(param[2] - 2) * (1 - u)^(2 * param[1] - 2) * 
      prod(param) - (-log(-(1 - u)^param[1] + 1))^(param[2] - 2) * (1 - u)^(2 * param[1] - 2) * param[1] - 
      (-log(-(1 - u)^param[1] + 1))^(param[2] - 1) * (1 - u)^(param[1] - 2) * param[1] - (-log(-(1 - 
      u)^param[1] + 1))^(param[2] - 1) * (1 - u)^(2 * param[1] - 2) + (-log(-(1 - u)^param[1] + 1))^(param[2] - 
      1) * (1 - u)^(param[1] - 2)))/((1 - u)^param[1] - 1)^2
  } else if (copula == 9) {
    out <- prod(param) * (1 - u)^(param[1] - 2) * (1 - (1 - u)^param[1])^(-2 - param[2]) * ((1 + prod(param)) * 
      (1 - u)^param[1] + param[1] - 1)
  } else if (copula == 10) {
    out <- (param[2]^2 * param[1] * ((1 - u * param[2])^(param[1] - 2) * param[1] + (1 - u * param[2])^(2 * 
      param[1] - 2) - (1 - u * param[2])^(param[1] - 2)))/(((1 - u * param[2])^param[1] - 1)^2)
  }
  
  return(out)
}


#################################################################
# cop.cdf  						#
# Input:							#
# u1,u2		data vectors					#
# param		copula parameters				#
# copula	copula family (7=BB1, 8=BB6, 9=BB7, 10=BB8)	#
# Output:							#
# out		copula						#
#################################################################

cop.cdf <- function(u1, u2, param, copula) {
  return(gen.inv(u = gen(u = u1, param = param, copula = copula) + gen(u = u2, param = param, copula = copula), 
    param = param, copula = copula))
}


#########################
# pdf for BB6  	#
#########################

bb6pdf <- function(u, v, th, de) {
  t1 <- 1 - u
  t2 <- t1^th
  t3 <- 1 - t2
  t4 <- log(t3)
  t5 <- (-t4)^de
  t12 <- 1/de
  t16 <- 1/th
  t32 <- de - 1
  t38 <- 2 * de
  t39 <- -1 + t38
  t40 <- (-t4)^t39
  t47 <- 3 * de - 1
  t50 <- (-t4)^t32
  t61 <- (-t4)^t47
  t90 <- (-t4)^t38
  # above depend on u and parameters only loop below for solving for conditional quantile
  t6 <- 1 - v
  t7 <- t6^th
  t8 <- 1 - t7
  t9 <- log(t8)
  t10 <- (-t9)^de
  t11 <- t5 + t10
  t13 <- t11^t12
  t14 <- exp(-t13)
  t15 <- 1 - t14
  t17 <- t15^t16
  t35 <- t11^(-2 * t32 * t12)
  t36 <- t35 * th
  t37 <- exp(t13)
  t42 <- (-t9)^t39
  t48 <- (-t9)^t47
  t53 <- t13 * de
  t56 <- (-t9)^t32
  t57 <- t37 * t50 * t56
  t59 <- t13 * th
  t78 <- t37 - 1
  t80 <- (t78 * t14)^t16
  t87 <- t78 * t78
  t93 <- (-t9)^t38
  # c21 = -t17*t13*t5*t2/t1/t3/t4/t11*t14/t15;
  pdf <- (2 * t36 * t37 * t40 * t42 + t36 * t37 * t48 * t50 + t53 * th * t57 - t59 * t57 + t36 * t37 * 
    t61 * t56 - 2 * t35 * t40 * t42 - t35 * t61 * t56 - t53 * th * t50 * t56 + t59 * t50 * t56 - t35 * 
    t48 * t50) * t80 * t7 * t2/t3/t8/t87/(t90 + 2 * t5 * t10 + t93)/t1/t6
  pdf
}


#################################################################################
# cop.pdf  								#
# Input:									#
# u1,u2		data vectors							#
# param		copula parameters						#
# copula	copula family (1,2,3,...14)					#
# Output:									#
# out		copula density							#
#################################################################################

cop.pdf <- function(u1, u2, param, copula) {
  if (copula == 7 | copula == 9 | copula == 10) {
    return(-gen.drv2(u = cop.cdf(u1 = u1, u2 = u2, param = param, copula = copula), param = param, copula = copula) * 
      gen.drv(u = u1, param = param, copula = copula) * gen.drv(u = u2, param = param, copula = copula)/gen.drv(u = cop.cdf(u1 = u1, 
      u2 = u2, param = param, copula = copula), param = param, copula = copula)^3)
  } else if (copula == 8) {
    return(bb6pdf(u1, u2, param[1], param[2]))
  } else if (copula == 17 | copula == 19 | copula == 20) {
    d1 <- 1 - u1
    d2 <- 1 - u2
    return(-gen.drv2(u = cop.cdf(u1 = d1, u2 = d2, param = param, copula = copula - 10), param = param, 
      copula = copula - 10) * gen.drv(u = d1, param = param, copula = copula - 10) * gen.drv(u = d2, 
      param = param, copula = copula - 10)/gen.drv(u = cop.cdf(u1 = d1, u2 = d2, param = param, copula = copula - 
      10), param = param, copula = copula - 10)^3)
  } else if (copula == 18) {
    d1 <- 1 - u1
    d2 <- 1 - u2
    return(bb6pdf(d1, d2, param[1], param[2]))
  } else if (copula == 27 | copula == 29 | copula == 30) {
    d1 <- 1 - u1
    d2 <- u2
    param <- -param
    return(-gen.drv2(u = cop.cdf(u1 = d1, u2 = d2, param = param, copula = copula - 20), param = param, 
      copula = copula - 20) * gen.drv(u = d1, param = param, copula = copula - 20) * gen.drv(u = d2, 
      param = param, copula = copula - 20)/gen.drv(u = cop.cdf(u1 = d1, u2 = d2, param = param, copula = copula - 
      20), param = param, copula = copula - 20)^3)
  } else if (copula == 28) {
    d1 <- 1 - u1
    d2 <- u2
    return(bb6pdf(d1, d2, -param[1], -param[2]))
  } else if (copula == 37 | copula == 39 | copula == 40) {
    d1 <- u1
    d2 <- 1 - u2
    param <- -param
    return(-gen.drv2(u = cop.cdf(u1 = d1, u2 = d2, param = param, copula = copula - 30), param = param, 
      copula = copula - 30) * gen.drv(u = d1, param = param, copula = copula - 30) * gen.drv(u = d2, 
      param = param, copula = copula - 30)/gen.drv(u = cop.cdf(u1 = d1, u2 = d2, param = param, copula = copula - 
      30), param = param, copula = copula - 30)^3)
  } else if (copula == 38) {
    d1 <- u1
    d2 <- 1 - u2
    return(bb6pdf(d1, d2, -param[1], -param[2]))
  } else if (copula == 0) {
    # independent
    return(rep(1, length(u1)))
  } else if (copula == 1) {
    # Gaussian
    t1 <- qnorm(p = u1)
    t2 <- qnorm(p = u2)
    rho <- param
    return(1/sqrt(1 - rho^2) * exp(-(rho^2 * (t1^2 + t2^2) - 2 * rho * t1 * t2)/(2 * (1 - rho^2))))
  } else if (copula == 2) {
    # t-copula
    rho <- param[1]
    nu <- param[2]
    
    t1 <- qt(u1, nu)
    t2 <- qt(u2, nu)
    return(1/(2 * pi * sqrt(1 - rho^2) * dt(t1, nu) * dt(t2, nu)) * (1 + (t1^2 + t2^2 - 2 * rho * t1 * 
      t2)/(nu * (1 - rho^2)))^(-(nu + 2)/2))
  } else if (copula == 3) {
    # Clayton
    theta <- param
    return((1 + theta) * (u1 * u2)^(-1 - theta) * (u1^(-theta) + u2^(-theta) - 1)^(-2 - 1/theta))
  } else if (copula == 4) {
    # Gumbel
    theta <- param
    t1 <- (-log(u1))^(theta) + (-log(u2))^(theta)
    t2 <- exp(-t1^(1/theta))
    return(t2/(u1 * u2) * t1^(-2 + 2/theta) * (log(u1) * log(u2))^(theta - 1) * (1 + (theta - 1) * t1^(-1/theta)))
  } else if (copula == 5) {
    # Frank
    theta <- param
    return((theta * (exp(theta) - 1) * exp(theta * u2 + theta * u1 + theta))/(exp(theta * u2 + theta * 
      u1) - exp(theta * u2 + theta) - exp(theta * u1 + theta) + exp(theta))^2)
  } else if (copula == 6) {
    # Joe
    theta <- param
    return(((1 - u1)^(theta) + (1 - u2)^(theta) - (1 - u1)^(theta) * (1 - u2)^(theta))^(1/(theta) - 2) * 
      (1 - u1)^(theta - 1) * (1 - u2)^(theta - 1) * (theta - 1 + (1 - u1)^(theta) + (1 - u2)^(theta) - 
      (1 - u1)^(theta) * (1 - u2)^(theta)))
  } else if (copula == 13) {
    # rotated Clayton (180)
    theta <- param
    d1 <- 1 - u1
    d2 <- 1 - u2
    return((1 + theta) * (d1 * d2)^(-1 - theta) * (d1^(-theta) + d2^(-theta) - 1)^(-2 - 1/theta))
  } else if (copula == 14) {
    # rotated Gumbel (180)
    theta <- param
    d1 <- 1 - u1
    d2 <- 1 - u2
    t1 <- (-log(d1))^(theta) + (-log(d2))^(theta)
    t2 <- exp(-t1^(1/theta))
    return(t2/(d1 * d2) * t1^(-2 + 2/theta) * (log(d1) * log(d2))^(theta - 1) * (1 + (theta - 1) * t1^(-1/theta)))
  } else if (copula == 16) {
    # rotated Joe (180)
    theta <- param
    d1 <- 1 - u1
    d2 <- 1 - u2
    return(((1 - d1)^(theta) + (1 - d2)^(theta) - (1 - d1)^(theta) * (1 - d2)^(theta))^(1/(theta) - 2) * 
      (1 - d1)^(theta - 1) * (1 - d2)^(theta - 1) * (theta - 1 + (1 - d1)^(theta) + (1 - d2)^(theta) - 
      (1 - d1)^(theta) * (1 - d2)^(theta)))
  } else if (copula == 23) {
    # rotated Clayton (90)
    theta <- -param
    d1 <- 1 - u1
    d2 <- u2
    return((1 + theta) * (d1 * d2)^(-1 - theta) * (d1^(-theta) + d2^(-theta) - 1)^(-2 - 1/theta))
  } else if (copula == 24) {
    # rotated Gumbel (90)
    theta <- -param
    d1 <- 1 - u1
    d2 <- u2
    t1 <- (-log(d1))^(theta) + (-log(d2))^(theta)
    t2 <- exp(-t1^(1/theta))
    return(t2/(d1 * d2) * t1^(-2 + 2/theta) * (log(d1) * log(d2))^(theta - 1) * (1 + (theta - 1) * t1^(-1/theta)))
  } else if (copula == 26) {
    # rotated Joe (90)
    theta <- -param
    d1 <- 1 - u1
    d2 <- u2
    return(((1 - d1)^(theta) + (1 - d2)^(theta) - (1 - d1)^(theta) * (1 - d2)^(theta))^(1/(theta) - 2) * 
      (1 - d1)^(theta - 1) * (1 - d2)^(theta - 1) * (theta - 1 + (1 - d1)^(theta) + (1 - d2)^(theta) - 
      (1 - d1)^(theta) * (1 - d2)^(theta)))
  } else if (copula == 33) {
    # rotaed Clayton (270)
    theta <- -param
    d1 <- u1
    d2 <- 1 - u2
    return((1 + theta) * (d1 * d2)^(-1 - theta) * (d1^(-theta) + d2^(-theta) - 1)^(-2 - 1/theta))
  } else if (copula == 34) {
    # rotated Gumbel (270)
    theta <- -param
    d1 <- u1
    d2 <- 1 - u2
    t1 <- (-log(d1))^(theta) + (-log(d2))^(theta)
    t2 <- exp(-t1^(1/theta))
    return(t2/(d1 * d2) * t1^(-2 + 2/theta) * (log(d1) * log(d2))^(theta - 1) * (1 + (theta - 1) * t1^(-1/theta)))
  } else if (copula == 36) {
    # rotated Joe (270)
    theta <- -param
    d1 <- u1
    d2 <- 1 - u2
    return(((1 - d1)^(theta) + (1 - d2)^(theta) - (1 - d1)^(theta) * (1 - d2)^(theta))^(1/(theta) - 2) * 
      (1 - d1)^(theta - 1) * (1 - d2)^(theta - 1) * (theta - 1 + (1 - d1)^(theta) + (1 - d2)^(theta) - 
      (1 - d1)^(theta) * (1 - d2)^(theta)))
  }
}



#########################################################
# Density of a meta-distribution with normal margins  #
# Input:						#
# x1,x2		vectors					#
# param		Copula parameter(s)			#
# copula	copula family				#
# Output:						#
# density						#
#########################################################

meta.dens <- function(x1, x2, param, copula, margins, margins.par) {
  if (margins == "norm") 
    return(cop.pdf(u1 = pnorm(x1), u2 = pnorm(x2), param = param, copula = copula) * dnorm(x1) * dnorm(x2)) else if (margins == "t") 
    return(cop.pdf(u1 = pt(x1, df = margins.par), u2 = pt(x2, df = margins.par), param = param, copula = copula) * 
      dt(x1, df = margins.par) * dt(x2, df = margins.par)) else if (margins == "unif") 
    return(cop.pdf(u1 = x1, u2 = x2, param = param, copula = copula)) else if (margins == "gamma") 
    return(cop.pdf(u1 = pgamma(x1, shape = margins.par[1], scale = margins.par[2]), u2 = pgamma(x2, shape = margins.par[1], 
      scale = margins.par[2]), param = param, copula = copula) * dgamma(x1, shape = margins.par[1], 
      scale = margins.par[2]) * dgamma(x2, shape = margins.par[1], scale = margins.par[2])) else if (margins == "exp") 
    return(cop.pdf(u1 = pexp(x1, rate = margins.par), u2 = pexp(x2, rate = margins.par), param = param, 
      copula = copula) * dexp(x1, rate = margins.par) * dexp(x2, rate = margins.par))
}



#################################################
# CopulaContour2D  			#
# Input:					#
# u1, u2	data-vectors			#
# bw		bandwidth			#
# size		number of grid point		#
# levels	Vector of contour levels	#
# family	copula family			#
# par		Copula parameter(s)		#
# Output:					#
# theo. and emp. contourplot			#
#################################################



BiCopMetaContour <- function(u1 = NULL, u2 = NULL, bw = 1, size = 100, levels = c(0.01, 0.05, 0.1, 0.15, 
  0.2), family = "emp", par = 0, par2 = 0, PLOT = TRUE, margins = "norm", margins.par = 0, xylim = NA, 
  ...) {
  if ((is.null(u1) == TRUE || is.null(u2) == TRUE) && family == "emp") 
    stop("'u1' and/or 'u2' not set or of length zero.")
  if (is.null(u1) != TRUE && (any(u1 > 1) || any(u1 < 0))) 
    stop("Data has be in the interval [0,1].")
  if (is.null(u2) != TRUE && (any(u2 > 1) || any(u2 < 0))) 
    stop("Data has be in the interval [0,1].")
  # if(length(u1)!=length(u2)) stop('Lengths of 'u1' and 'u2' do not match.')
  if (!(family %in% c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 13, 14, 16, 17, 18, 19, 20, 23, 24, 26, 27, 28, 
    29, 30, 33, 34, 36, 37, 38, 39, 40, "emp"))) 
    stop("Copula family not implemented.")
  if (family %in% c(2, 7, 8, 9, 10, 17, 18, 19, 20, 27, 28, 29, 30, 37, 38, 39, 40) && par2 == 0) 
    stop("For t-, BB1 and BB7 copulas, 'par2' must be set.")
  if (family %in% c(1, 3, 4, 5, 6, 13, 14, 16, 23, 24, 26, 33, 34, 36) && length(par) < 1) 
    stop("'par' not set.")
  
  # size sollte nicht zu gross sein
  if (size > 1000) 
    stop("Size parameter should not be greater than 1000. Otherwise computational time and memory space are too large.")
  if (size < 50) 
    stop("Size parameter should not be smaller than 50.")
  
  # bw richtig
  if (bw < 1) 
    stop("The bandwidth parameter 'bw' should be greater or equal to 1.")
  if (bw > 5) 
    stop("The bandwidth parameter 'bw' should not be greater than 5.")
  
  # Parameterbereiche abfragen
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
  
  
  if (PLOT != TRUE && PLOT != FALSE) 
    stop("The parameter 'PLOT' has to be set to 'TRUE' or 'FALSE'.")
  
  if (margins != "norm" && margins != "t" && margins != "exp" && margins != "gamma" && margins != "unif") 
    stop("The function only supports Gaussian ('norm'), Student t ('t'), exponential ('exp'), Gamma ('gamma') and uniform ('unif') margins.")
  
  if (margins == "t" && margins.par <= 0) 
    stop("The degrees of freedom parameter for the Student t margins has to positive.")
  if (margins == "Gamma" && length(margins.par) != 2) 
    stop("For Gamma margins two parameters are required in 'margins.par'.")
  if (margins == "exp" && margins.par == 0) 
    stop("Exponential margins require one parameter in 'margins.par'.")
  if (margins == "unif" && family == 0) 
    stop("The combination independence copula and uniform margins is not possible because all z-values are equal.")
  
  if (is.null(u1) && is.null(u2) && family != "emp") {
    # theoretischer contourplot
    u1 <- runif(1000)
    u2 <- runif(1000)
  }
  
  if (!is.na(xylim) && length(xylim) != 2) 
    stop("'xylim' has to be a vector of length 2.")
  
  if (margins == "norm") {
    x1 <- qnorm(p = u1)
    x2 <- qnorm(p = u2)
    if (is.na(xylim)) 
      xylim <- c(-3, 3)
  } else if (margins == "t") {
    x1 <- qt(p = u1, df = margins.par)
    x2 <- qt(p = u2, df = margins.par)
    if (is.na(xylim)) 
      xylim <- c(-3, 3)
  } else if (margins == "exp") {
    x1 <- qexp(p = u1, rate = margins.par)
    x2 <- qexp(p = u2, rate = margins.par)
    if (is.na(xylim)) 
      xylim <- c(0, 5)
  } else if (margins == "gamma") {
    x1 <- qgamma(p = u1, shape = margins.par[1], scale = margins.par[2])
    x2 <- qgamma(p = u2, shape = margins.par[1], scale = margins.par[2])
    if (is.na(xylim)) 
      xylim <- c(0, 5)
  } else if (margins == "unif") {
    x1 <- u1
    x2 <- u2
    if (is.na(xylim)) 
      xylim <- c(0, 1)
  }
  
  x <- y <- seq(from = xylim[1], to = xylim[2], length.out = size)
  
  if (family != "emp") {
    if (family %in% c(2, 7, 8, 9, 10, 17, 18, 19, 20, 27, 28, 29, 30, 37, 38, 39, 40)) 
      z <- matrix(data = meta.dens(x1 = rep(x = x, each = size), x2 = rep(x = y, times = size), param = c(par, 
        par2), copula = family, margins = margins, margins.par = margins.par), nrow = size, byrow = TRUE) else z <- matrix(data = meta.dens(x1 = rep(x = x, each = size), x2 = rep(x = y, times = size), param = par, 
      copula = family, margins = margins, margins.par = margins.par), nrow = size, byrow = TRUE)
  } else {
    # empirical
    bw1 <- bw * bandwidth.nrd(x1)
    bw2 <- bw * bandwidth.nrd(x2)
    
    kd.est <- kde2d(x = x1, y = x2, h = c(bw1, bw2), n = size)
    
    x <- kd.est$x
    y <- kd.est$y
    z <- kd.est$z
  }
  
  if (PLOT) {
    contour(x = x, y = y, z = z, levels = levels, ylim = xylim, xlim = xylim, ...)
  } else {
    out <- list()
    out$x <- x
    out$y <- y
    out$z <- z
    
    return(out)
  }
} 
