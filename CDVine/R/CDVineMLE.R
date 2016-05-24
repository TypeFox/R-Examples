CDVineMLE <- function(data, family, start = NULL, start2 = NULL, type, maxit = 200, max.df = 30,
                      max.BB = list(BB1 = c(5, 6), BB6 = c(6, 6), BB7 = c(5, 6), BB8 = c(6, 1)), ...) {
  if (type == "CVine") 
    type <- 1 else if (type == "DVine") 
    type <- 2
  if (type != 1 & type != 2) 
    stop("Vine model not implemented.")
  
  d <- dim(data)[2]
  T <- dim(data)[1]
  Maxiter <- floor(maxit)
  if (any(data > 1) || any(data < 0)) 
    stop("Data has be in the interval [0,1].")
  
  if (max.df <= 2) 
    stop("The upper bound for the degrees of freedom parameter has to be larger than 2.")
  if (!is.list(max.BB)) 
    stop("'max.BB' has to be a list.")
  if (max.BB$BB1[1] < 0.001) 
    stop("The upper bound for the first parameter of the BB1 copula should be greater than 0.001 (lower bound for estimation).")
  if (max.BB$BB1[2] < 1.001) 
    stop("The upper bound for the second parameter of the BB1 copula should be greater than 1.001 (lower bound for estimation).")
  if (max.BB$BB6[1] < 1.001) 
    stop("The upper bound for the first parameter of the BB6 copula should be greater than 1.001 (lower bound for estimation).")
  if (max.BB$BB6[2] < 1.001) 
    stop("The upper bound for the second parameter of the BB6 copula should be greater than 1.001 (lower bound for estimation).")
  if (max.BB$BB7[1] < 1.001) 
    stop("The upper bound for the first parameter of the BB7 copula should be greater than 1.001 (lower bound for estimation).")
  if (max.BB$BB7[2] < 0.001) 
    stop("The upper bound for the second parameter of the BB7 copula should be greater than 0.001 (lower bound for estimation).")
  if (max.BB$BB8[1] < 1.001) 
    stop("The upper bound for the first parameter of the BB1 copula should be greater than 0.001 (lower bound for estimation).")
  if (max.BB$BB8[2] < 0.001 || max.BB$BB8[2] > 1) 
    stop("The upper bound for the second parameter of the BB1 copula should be in the interval [0,1].")
  
  if (length(family) != d * (d - 1)/2) 
    stop("Number of copula families incorrect.")
  
  if (!is.null(start) && length(start) != (d * (d - 1)/2)) 
    stop("Length of the vector 'start' is incorrect.")
  if (!is.null(start2) && length(start2) != (d * (d - 1)/2)) 
    stop("Length of the vector 'start2' is incorrect.")
  if ((is.null(start) && !is.null(start2)) || (!is.null(start) && is.null(start2))) 
    stop("If one of the starting parameter vectors is set, the other one has to be set, too.")
  # Sicherheitsabfrage
  for (i in 1:(d * (d - 1)/2)) {
    if (!(family[i] %in% c(0, 1:10, 13, 14, 16:20, 23, 24, 26:30, 33, 34, 36:40))) 
      stop("Copula family not implemented.")
    if (!is.null(start)) {
      # Parameterbereiche abfragen
      if ((family[i] == 1 || family[i] == 2) && abs(start[i]) >= 1) 
        stop("The parameter of the Gaussian and t-copula has to be in the interval (-1,1).")
      if (family[i] == 2 && start2[i] <= 2) 
        stop("The degrees of freedom parameter of the t-copula has to be larger than 2.")
      if ((family[i] == 3 || family[i] == 13) && start[i] <= 0) 
        stop("The parameter of the Clayton copula has to be positive.")
      if ((family[i] == 4 || family[i] == 14) && start[i] < 1) 
        stop("The parameter of the Gumbel copula has to be in the interval [1,oo).")
      if ((family[i] == 6 || family[i] == 16) && start[i] <= 1) 
        stop("The copula parameter of the Joe copula has to be in the interval (1,oo).")
      if (family[i] == 5 && start[i] == 0) 
        stop("The parameter of the Frank copula has to be unequal to 0.")
      if ((family[i] == 7 || family[i] == 17) && start[i] <= 0) 
        stop("The first parameter of the BB1 copula has to be positive.")
      if ((family[i] == 7 || family[i] == 17) && start2[i] < 1) 
        stop("The second parameter of the BB1 copula has to be in the interval [1,oo).")
      if ((family[i] == 8 || family[i] == 18) && start[i] <= 0) 
        stop("The first parameter of the BB6 copula has to be in the interval [1,oo).")
      if ((family[i] == 8 || family[i] == 18) && start2[i] < 1) 
        stop("The second parameter of the BB6 copula has to be in the interval [1,oo).")
      if ((family[i] == 9 || family[i] == 19) && start[i] < 1) 
        stop("The first parameter of the BB7 copula has to be in the interval [1,oo).")
      if ((family[i] == 9 || family[i] == 19) && start2[i] <= 0) 
        stop("The second parameter of the BB7 copula has to be positive.")
      if ((family[i] == 10 || family[i] == 20) && start[i] < 1) 
        stop("The first parameter of the BB8 copula has to be in the interval [1,oo).")
      if ((family[i] == 10 || family[i] == 20) && (start2[i] <= 0 || start2[i] > 1)) 
        stop("The second parameter of the BB8 copula has to be in the interval (0,1].")
      if ((family[i] == 23 || family[i] == 33) && start[i] >= 0) 
        stop("The parameter of the rotated Clayton copula has to be negative.")
      if ((family[i] == 24 || family[i] == 34) && start[i] > -1) 
        stop("The parameter of the rotated Gumbel copula has to be in the interval (-oo,-1].")
      if ((family[i] == 26 || family[i] == 36) && start[i] >= -1) 
        stop("The parameter of the rotated Joe copula has to be in the interval (-oo,-1).")
      if ((family[i] == 27 || family[i] == 37) && start[i] >= 0) 
        stop("The first parameter of the rotated BB1 copula has to be negative.")
      if ((family[i] == 27 || family[i] == 37) && start2[i] > -1) 
        stop("The second parameter of the rotated BB1 copula has to be in the interval (-oo,-1].")
      if ((family[i] == 28 || family[i] == 38) && start[i] >= 0) 
        stop("The first parameter of the rotated BB6 copula has to be in the interval (-oo,-1].")
      if ((family[i] == 28 || family[i] == 38) && start2[i] > -1) 
        stop("The second parameter of the rotated BB6 copula has to be in the interval (-oo,-1].")
      if ((family[i] == 29 || family[i] == 39) && start[i] > -1) 
        stop("The first parameter of the rotated BB7 copula has to be in the interval (-oo,-1].")
      if ((family[i] == 29 || family[i] == 39) && start2[i] >= 0) 
        stop("The second parameter of the rotated BB7 copula has to be negative.")
      if ((family[i] == 30 || family[i] == 40) && start[i] > -1) 
        stop("The first parameter of the rotated BB8 copula has to be in the interval (-oo,-1].")
      if ((family[i] == 30 || family[i] == 40) && (start2[i] >= 0 || start2[i] < (-1))) 
        stop("The second parameter of the rotated BB8 copula has to be in the interval [-1,0).")
    }
  }
  if (T < 2) 
    stop("Number of observations has to be at least 2.")
  if (d < 2) 
    stop("Dimension has to be at least 2.")
  if (Maxiter < 1) 
    stop("'maxit' has to be greater than zero.")
  
  dd <- d * (d - 1)/2
  tt <- sum(family == 2) + sum(family == 7) + sum(family == 8) + sum(family == 9) + sum(family == 10) + 
    sum(family == 17) + sum(family == 18) + sum(family == 19) + sum(family == 20) + sum(family == 27) + 
    sum(family == 28) + sum(family == 29) + sum(family == 30) + sum(family == 37) + sum(family == 38) + 
    sum(family == 39) + sum(family == 40)
  start_par <- list()
  if (is.null(start)) 
    start_par <- CDVineSeqEst(data, family, type, max.df = max.df) else {
    start_par$par <- start
    start_par$par2 <- start2
  }
  
  
  
  parm0 <- start_par$par[(family != 0)]
  for (k in 1:dd) {
    if (family[k] %in% c(2, 7, 8, 9, 10, 17, 18, 19, 20, 27, 28, 29, 30, 37, 38, 39, 40)) {
      parm0 <- c(parm0, start_par$par2[k])
    }
    
  }
  
  param1 <- parm0
  pscale1 <- numeric()
  for (i in 1:dd) {
    pscale1[i] <- ifelse(family[i] == 1, 0.01, 1)
  }
  pscale <- c(pscale1[family != 0], rep(1, tt))
  
  
  func <- function(parm, data, family, type) {
    for (k in 1:dd) {
      # Die independent copula wieder aufnehmen
      if (family[k] == 0) {
        if (k == 1) {
          parm <- c(0, parm)
        } else if (k > length(parm)) {
          parm <- c(parm, 0)
        } else {
          parm <- c(parm[1:(k - 1)], 0, parm[k:length(parm)])
        }
      }
    }
    
    
    parm1 <- parm[1:dd]
    nu <- numeric()
    kk <- 1
    for (k in 1:dd) {
      if (family[k] %in% c(2, 7, 8, 9, 10, 17, 18, 19, 20, 27, 28, 29, 30, 37, 38, 39, 40)) {
        nu[k] <- parm[dd + kk]
        kk <- kk + 1
      } else {
        nu[k] <- 0
      }
    }
    CDVineLogLik(data, family, parm1, nu, type)$loglik
  }
  
  # Grenzen setzen fuer die Optimierung
  l <- rep(0, dd)
  u <- rep(1, dd)
  for (j in 1:dd) {
    if (family[j] == 1 | family[j] == 2) {
      l[j] <- -0.9999
      u[j] <- 0.9999
    } else if (family[j] == 3 | family[j] == 13) {
      l[j] <- 1e-04
      u[j] <- Inf
    } else if (family[j] == 4 | family[j] == 14) {
      l[j] <- 1.0001
      u[j] <- Inf
    } else if (family[j] == 5) {
      l[j] <- -Inf
      u[j] <- Inf
    } else if (family[j] == 6 | family[j] == 16) {
      l[j] <- 1.0001
      u[j] <- Inf
    } else if (family[j] == 7 | family[j] == 17) {
      l[j] <- 0.001
      u[j] <- max.BB$BB1[1]
    } else if (family[j] == 8 | family[j] == 18) {
      l[j] <- 1.001
      u[j] <- max.BB$BB6[1]
    } else if (family[j] == 9 | family[j] == 19) {
      l[j] <- 1.001
      u[j] <- max.BB$BB7[1]
    } else if (family[j] == 10 | family[j] == 20) {
      l[j] <- 1.001
      u[j] <- max.BB$BB8[1]
    } else if (family[j] == 23 | family[j] == 33) {
      l[j] <- -Inf
      u[j] <- -1e-04
    } else if (family[j] == 24 | family[j] == 34) {
      l[j] <- -Inf
      u[j] <- -1.0001
    } else if (family[j] == 26 | family[j] == 36) {
      l[j] <- -Inf
      u[j] <- -1.0001
    } else if (family[j] == 27 | family[j] == 37) {
      l[j] <- -max.BB$BB1[1]
      u[j] <- -0.001
    } else if (family[j] == 28 | family[j] == 38) {
      l[j] <- -max.BB$BB6[1]
      u[j] <- -1.001
    } else if (family[j] == 29 | family[j] == 39) {
      l[j] <- -max.BB$BB7[1]
      u[j] <- -1.001
    } else if (family[j] == 30 | family[j] == 40) {
      l[j] <- -max.BB$BB8[1]
      u[j] <- -1.001
    } else if (family[j] == 0) {
      # independence copula => no parameter
      l[j] <- 0
      u[j] <- 0
    } else stop("Copula family not implemented.")
  }
  
  l <- l[family != 0]
  u <- u[family != 0]
  
  for (j in 1:dd) {
    if (family[j] == 2) {
      l <- c(l, 2.0001)
      u <- c(u, max.df)
    } else if (family[j] == 7 | family[j] == 17) {
      l <- c(l, 1.001)
      u <- c(u, max.BB$BB1[2])
    } else if (family[j] == 8 | family[j] == 18) {
      l <- c(l, 1.001)
      u <- c(u, max.BB$BB6[2])
    } else if (family[j] == 9 | family[j] == 19) {
      l <- c(l, 0.001)
      u <- c(u, max.BB$BB7[2])
    } else if (family[j] == 10 | family[j] == 20) {
      l <- c(l, 0.001)
      u <- c(u, max.BB$BB8[2])
    } else if (family[j] == 27 | family[j] == 37) {
      l <- c(l, -max.BB$BB1[2])
      u <- c(u, -1.001)
    } else if (family[j] == 28 | family[j] == 38) {
      l <- c(l, -max.BB$BB6[2])
      u <- c(u, -1.001)
    } else if (family[j] == 29 | family[j] == 39) {
      l <- c(l, -max.BB$BB7[2])
      u <- c(u, -0.001)
    } else if (family[j] == 30 | family[j] == 40) {
      l <- c(l, -max.BB$BB8[2])
      u <- c(u, -0.001)
    }
  }
  
  if (!exists("factr")) 
    factr <- 1e+08
  
  
  out <- optim(par = param1, fn = func, gr = NULL, data = data, family = family, type = type, method = "L-BFGS-B", 
    lower = l, upper = u, control = list(fnscale = -1, maxit = Maxiter, parscale = pscale, factr = factr, 
      ...))
  ## fnscale=-1 turns the problem into a maximization problem, maxit defines the maximal number of
  ## iterations
  
  for (kk in 1:dd) {
    # Die independent copula wieder einfuegen
    if (family[kk] == 0) {
      if (kk == 1) {
        out$par <- c(0, out$par)
      } else if (kk > length(out$par)) {
        out$par <- c(out$par, 0)
      } else {
        out$par <- c(out$par[1:(kk - 1)], 0, out$par[kk:length(out$par)])
      }
    }
  }
  
  
  
  # Parameter ordentlich ausgeben
  theta <- out$par[1:dd]
  nu <- numeric()
  kk <- 1
  for (k in 1:dd) {
    if (family[k] %in% c(2, 7, 8, 9, 10, 17, 18, 19, 20, 27, 28, 29, 30, 37, 38, 39, 40)) {
      nu[k] <- out$par[dd + kk]
      kk <- kk + 1
    } else {
      nu[k] <- 0
    }
  }
  
  out2 <- list()
  out2$par <- theta
  out2$par2 <- nu
  out2$loglik <- out$value
  out2$counts <- out$counts
  out2$convergence <- out$convergence
  out2$message <- out$message
  
  return(out2)
}






# Parametertransformationen

trafo1 <- function(p, fam, par) {
  if (fam == 1 || (fam == 2 && par == 1)) {
    pNew <- atanh(p)
  } else if (fam == 0 || fam == 5) {
    pNew <- p
  } else if (fam == 3 || (fam == 7 && par == 1) || (fam == 9 && par == 2) || fam == 8 || fam == 13) {
    pNew <- log(p)
  } else if (fam == 4 || fam == 14 || (fam == 2 && par == 2) || fam == 6 || (fam == 7 && par == 2) || (fam == 
    9 && par == 1) || fam == 16) {
    pNew <- log(p - 1)
  } else if (fam == 23 || fam == 33) {
    pNew <- log(-p)
  } else if (fam == 24 || fam == 26 || fam == 34 || fam == 36) {
    pNew <- log(-p + 1)
  }
  return(pNew)
}

trafo2 <- function(p, fam, par) {
  if (fam == 1 || (fam == 2 && par == 1)) {
    pNew <- tanh(p)
  } else if (fam == 0 || fam == 5) {
    pNew <- p
  } else if (fam == 3 || (fam == 7 && par == 1) || (fam == 9 && par == 2) || fam == 8 || fam == 13) {
    pNew <- exp(p)
    if (pNew %in% c(NA, NaN, Inf)) 
      pNew <- 1e+308 else if (pNew == -Inf) 
      pNew <- -1e+308
    if (fam == 7 && pNew < 1e-04) 
      pNew <- 1e-04
  } else if ((fam == 2 && par == 2) || fam == 4 || fam == 14 || fam == 6 || (fam == 7 && par == 2) || (fam == 
    9 && par == 1) || fam == 16) {
    pNew <- exp(p) + 1
    if (pNew %in% c(NA, NaN, Inf)) 
      pNew <- 1e+308 else if (pNew == -Inf) 
      pNew <- -1e+308
  } else if (fam == 23 || fam == 33) {
    pNew <- -exp(p)
    if (pNew %in% c(NA, NaN, -Inf)) 
      pNew <- -1e+308 else if (pNew == Inf) 
      pNew <- 1e+308
  } else if (fam == 24 || fam == 26 || fam == 34 || fam == 36) {
    pNew <- -exp(p) - 1
    if (pNew %in% c(NA, NaN, -Inf)) 
      pNew <- -1e+308 else if (pNew == Inf) 
      pNew <- 1e+308
  }
  return(pNew)
}

trafo2der <- function(p, fam, par) {
  if (fam == 1 || (fam == 2 && par == 1)) {
    pNew <- 1 - tanh(p)^2
  } else if (fam == 0 || fam == 5) {
    pNew <- 1
  } else if (fam == 3 || (fam == 7 && par == 1) || (fam == 9 && par == 2) || fam == 8 || fam == 13) {
    pNew <- exp(p)
    if (pNew %in% c(NA, NaN, Inf)) 
      pNew <- 1e+308 else if (pNew == -Inf) 
      pNew <- -1e+308
  } else if ((fam == 2 && par == 2) || fam == 4 || fam == 14 || fam == 6 || (fam == 7 && par == 2) || (fam == 
    9 && par == 1) || fam == 16) {
    pNew <- exp(p)
    if (pNew %in% c(NA, NaN, Inf)) 
      pNew <- 1e+308 else if (pNew == -Inf) 
      pNew <- -1e+308
  } else if (fam == 23 || fam == 33) {
    pNew <- -exp(p)
    if (pNew %in% c(NA, NaN, -Inf)) 
      pNew <- -1e+308 else if (pNew == Inf) 
      pNew <- 1e+308
  } else if (fam == 24 || fam == 26 || fam == 34 || fam == 36) {
    pNew <- -exp(p)
    if (pNew %in% c(NA, NaN, -Inf)) 
      pNew <- -1e+308 else if (pNew == Inf) 
      pNew <- 1e+308
  }
  return(pNew)
} 
