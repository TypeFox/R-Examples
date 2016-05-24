CDVineSeqEst <- function(data, family, type, method = "mle", se = FALSE, max.df = 30,
                         max.BB = list(BB1 = c(5, 6), BB6 = c(6, 6), BB7 = c(5, 6), BB8 = c(6, 1)),
                         progress = FALSE) {
  # Function that estimates the parameter(s)
  # by stages, It can be used as starting values
  #---------------------------------------------------------
  # INPUT:
  #   data      Data for which to estimate parameter
  #   family  The array definig the copulas in the pcc copula construction
  #   type	Type of vine (1=canonical vine, 2=D-Vine)
  # OUTPUT:
  #   out     Estimated Parameters
  #----------------------------------------------------------
  # Author: Carlos Almeida and Ulf Schepsmeier
  # Date: 2009-11-11
  # Update: 2011-03-01 (by Ulf Schepsmeier)
  #---------------------------------------------------------------
  
  if (type == "CVine") 
    type <- 1 else if (type == "DVine") 
    type <- 2
  if (type != 1 & type != 2) 
    stop("Vine model not implemented.")
  
  if (method != "mle" && method != "itau") 
    stop("Estimation method has to be either 'mle' or 'itau'.")
  if (is.logical(se) == FALSE) 
    stop("'se' has to be a logical variable (TRUE or FALSE).")
  
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
  
  
  n <- dim(data)[1]
  d <- dim(data)[2]
  
  # Sicherheitsabfragen
  if (d < 3) 
    stop("Dimension has to be at least 3.")
  if (n < 2) 
    stop("Number of observations has to be at least 2.")
  if (any(data > 1) || any(data < 0)) 
    stop("Data has be in the interval [0,1].")
  if (length(family) != d * (d - 1)/2) 
    stop("Number of copula families incorrect.")
  for (i in 1:(d * (d - 1)/2)) {
    if (!(family[i] %in% c(0, 1:10, 13, 14, 16:20, 23, 24, 26:30, 33, 34, 36:40))) 
      stop("Copula family not implemented.")
  }
  
  nuMat <- matrix(0, nrow = d - 1, ncol = d - 1)
  rhoMat <- matrix(0, nrow = d - 1, ncol = d - 1)
  
  w <- matrix(0, d - 1, d - 1)
  k <- 1
  for (i in 1:(d - 1)) {
    for (j in 1:(d - i)) {
      w[i, j] <- family[k]
      k <- k + 1
    }
  }
  
  if (se == TRUE) {
    seMat1 <- matrix(0, nrow = d - 1, ncol = d - 1)
    seMat2 <- matrix(0, nrow = d - 1, ncol = d - 1)
  }
  
  if (type == 1) {
    # C-Vine
    v <- array(0, c(d - 1, d - 1, n))
    
    for (i in 1:(d - 1)) {
      if (w[1, i] %in% c(2, 7, 8, 9, 10, 17, 18, 19, 20, 27, 28, 29, 30, 37, 38, 39, 40)) {
        if (progress == TRUE) 
          message(1, ",", i + 1)
        par.out <- BiCopEst(data[, 1], data[, i + 1], w[1, i], method, se, max.df, max.BB = max.BB)
        # par1 <- par.out$par
        rhoMat[1, i] <- par.out$par
        nuMat[1, i] <- par.out$par2
        if (se == TRUE) {
          # se1 <- par.out$se
          seMat1[1, i] <- par.out$se
          seMat2[1, i] <- par.out$se2
        }
      } else {
        if (progress == TRUE) 
          message(1, ",", i + 1)
        par.out <- BiCopEst(data[, 1], data[, i + 1], w[1, i], method, se, max.df, max.BB = max.BB)
        rhoMat[1, i] <- par.out$par
        if (se == TRUE) {
          seMat1[1, i] <- par.out$se
        }
      }
      
      v[1, i, ] <- .C("Hfunc1", as.integer(w[1, i]), as.integer(n), as.double(data[, i + 1]), as.double(data[, 
        1]), as.double(rhoMat[1, i]), as.double(nuMat[1, i]), as.double(rep(0, n)), PACKAGE = "CDVine")[[7]]
    }
    for (j in 2:(d - 1)) {
      for (i in 1:(d - j)) {
        if (w[j, i] %in% c(2, 7, 8, 9, 10, 17, 18, 19, 20, 27, 28, 29, 30, 37, 38, 39, 40)) {
          if (progress == TRUE) 
          message(j, ",", j + i, "|", paste(1:(j - 1), collapse = ","))
          par.out <- BiCopEst(v[j - 1, 1, ], v[j - 1, i + 1, ], w[j, i], method, se, max.df, max.BB = max.BB)
          # par1 <- par.out$par
          rhoMat[j, i] <- par.out$par
          nuMat[j, i] <- par.out$par2
          if (se == TRUE) {
          # se1 <- par.out$se
          seMat1[j, i] <- par.out$se
          seMat2[j, i] <- par.out$se2
          }
        } else {
          if (progress == TRUE) 
          message(j, ",", j + i, "|", paste(1:(j - 1), collapse = ","))
          par.out <- BiCopEst(v[j - 1, 1, ], v[j - 1, i + 1, ], w[j, i], method, se, max.df, max.BB = max.BB)
          rhoMat[j, i] <- par.out$par
          if (se == TRUE) {
          seMat1[j, i] <- par.out$se
          }
        }
        
        if (j < (d - 1)) {
          v[j, i, ] <- .C("Hfunc1", as.integer(w[j, i]), as.integer(n), as.double(v[j - 1, i + 1, 
          ]), as.double(v[j - 1, 1, ]), as.double(rhoMat[j, i]), as.double(nuMat[j, i]), as.double(rep(0, 
          n)), PACKAGE = "CDVine")[[7]]
        }
      }
    }
  } else {
    # D-Vine
    v <- array(0, c(d, 2 * d - 4, n))
    
    for (i in 1:(d - 1)) {
      if (w[1, i] %in% c(2, 7, 8, 9, 10, 17, 18, 19, 20, 27, 28, 29, 30, 37, 38, 39, 40)) {
        if (progress == TRUE) 
          message(i, ",", i + 1)
        par.out <- BiCopEst(data[, i], data[, i + 1], w[1, i], method, se, max.df, max.BB = max.BB)
        # par1 <- par.out$par
        rhoMat[1, i] <- par.out$par
        nuMat[1, i] <- par.out$par2
        if (se == TRUE) {
          # se1 <- par.out$se
          seMat1[1, i] <- par.out$se
          seMat2[1, i] <- par.out$se2
        }
      } else {
        if (progress == TRUE) 
          message(i, ",", i + 1)
        par.out <- BiCopEst(data[, i], data[, i + 1], w[1, i], method, se, max.df, max.BB = max.BB)
        rhoMat[1, i] <- par.out$par
        if (se == TRUE) {
          seMat1[1, i] <- par.out$se
        }
      }
      
    }
    v[1, 1, ] <- .C("Hfunc2", as.integer(w[1, 1]), as.integer(n), as.double(data[, 1]), as.double(data[, 
      2]), as.double(rhoMat[1, 1]), as.double(nuMat[1, 1]), as.double(rep(0, n)), PACKAGE = "CDVine")[[7]]
    if (d > 3) {
      for (k in 1:(d - 3)) {
        v[1, 2 * k, ] <- .C("Hfunc1", as.integer(w[1, k + 1]), as.integer(n), as.double(data[, k + 
          2]), as.double(data[, k + 1]), as.double(rhoMat[1, k + 1]), as.double(nuMat[1, k + 1]), 
          as.double(rep(0, n)), PACKAGE = "CDVine")[[7]]
        v[1, 2 * k + 1, ] <- .C("Hfunc2", as.integer(w[1, k + 1]), as.integer(n), as.double(data[, 
          k + 1]), as.double(data[, k + 2]), as.double(rhoMat[1, k + 1]), as.double(nuMat[1, k + 
          1]), as.double(rep(0, n)), PACKAGE = "CDVine")[[7]]
      }
    }
    v[1, 2 * d - 4, ] <- .C("Hfunc1", as.integer(w[1, d - 1]), as.integer(n), as.double(data[, d]), as.double(data[, 
      d - 1]), as.double(rhoMat[1, d - 1]), as.double(nuMat[1, d - 1]), as.double(rep(0, n)), PACKAGE = "CDVine")[[7]]
    for (j in 2:(d - 1)) {
      for (i in 1:(d - j)) {
        if (w[j, i] %in% c(2, 7, 8, 9, 10, 17, 18, 19, 20, 27, 28, 29, 30, 37, 38, 39, 40)) {
          if (progress == TRUE) 
          message(i, ",", i + j, "|", paste((i + 1):(i + j - 1), collapse = ","))
          par.out <- BiCopEst(v[j - 1, 2 * i - 1, ], v[j - 1, 2 * i, ], w[j, i], method, se, max.df, 
          max.BB = max.BB)
          # par1 <- par.out$par
          rhoMat[j, i] <- par.out$par
          nuMat[j, i] <- par.out$par2
          if (se == TRUE) {
          # se1 <- par.out$se
          seMat1[j, i] <- par.out$se
          seMat2[j, i] <- par.out$se2
          }
        } else {
          if (progress == TRUE) 
          message(i, ",", i + j, "|", paste((i + 1):(i + j - 1), collapse = ","))
          par.out <- BiCopEst(v[j - 1, 2 * i - 1, ], v[j - 1, 2 * i, ], w[j, i], method, se, max.df, 
          max.BB = max.BB)
          rhoMat[j, i] <- par.out$par
          if (se == TRUE) {
          seMat1[j, i] <- par.out$se
          }
        }
        
      }
      v[j, 1, ] <- .C("Hfunc2", as.integer(w[j, 1]), as.integer(n), as.double(v[j - 1, 1, ]), as.double(v[j - 
        1, 2, ]), as.double(rhoMat[j, 1]), as.double(nuMat[j, 1]), as.double(rep(0, n)), PACKAGE = "CDVine")[[7]]
      if (d > 4 & (d - j - 2) > 0) {
        for (i in 1:(d - j - 2)) {
          v[j, 2 * i, ] <- .C("Hfunc1", as.integer(w[j, i + 1]), as.integer(n), as.double(v[j - 1, 
          2 * i + 2, ]), as.double(v[j - 1, 2 * i + 1, ]), as.double(rhoMat[j, i + 1]), as.double(nuMat[j, 
          i + 1]), as.double(rep(0, n)), PACKAGE = "CDVine")[[7]]
          v[j, 2 * i + 1, ] <- .C("Hfunc2", as.integer(w[j, i + 1]), as.integer(n), as.double(v[j - 
          1, 2 * i + 1, ]), as.double(v[j - 1, 2 * i + 2, ]), as.double(rhoMat[j, i + 1]), as.double(nuMat[j, 
          i + 1]), as.double(rep(0, n)), PACKAGE = "CDVine")[[7]]
        }
      }
      v[j, 2 * d - 2 * j - 2, ] <- .C("Hfunc1", as.integer(w[j, d - j]), as.integer(n), as.double(v[j - 
        1, 2 * d - 2 * j, ]), as.double(v[j - 1, 2 * d - 2 * j - 1, ]), as.double(rhoMat[j, d - j]), 
        as.double(nuMat[j, d - j]), as.double(rep(0, n)), PACKAGE = "CDVine")[[7]]
    }
  }
  
  theta0 <- rep(0, d * (d - 1)/2)
  se1 <- rep(0, d * (d - 1)/2)
  k <- 1
  for (j in 1:(d - 1)) {
    for (i in 1:(d - j)) {
      theta0[k] <- rhoMat[j, i]
      if (se == TRUE) 
        se1[k] <- seMat1[j, i]
      k <- k + 1
    }
  }
  tt <- sum(family == 2) + sum(family == 7) + sum(family == 8) + sum(family == 9)
  nu0 <- rep(0, tt)
  se0 <- rep(0, tt)
  k <- 1
  for (j in 1:(d - 1)) {
    for (i in 1:(d - j)) {
      if (w[j, i] %in% c(2, 7, 8, 9, 10, 17, 18, 19, 20, 27, 28, 29, 30, 37, 38, 39, 40)) {
        nu0[k] <- nuMat[j, i]
        if (se == TRUE) 
          se0[k] <- seMat2[j, i]
        k <- k + 1
      }
    }
  }
  
  nu1 <- numeric()
  se2 <- numeric()
  kk <- 1
  dd <- d * (d - 1)/2
  for (k in 1:dd) {
    if (family[k] %in% c(2, 7, 8, 9, 10, 17, 18, 19, 20, 27, 28, 29, 30, 37, 38, 39, 40)) {
      nu1[k] <- nu0[kk]
      if (se == TRUE) 
        se2[k] <- se0[kk]
      kk <- kk + 1
    } else {
      nu1[k] <- 0
      if (se == TRUE) 
        se2[k] <- 0
    }
  }
  
  if (se == TRUE) 
    out <- list(par = theta0, par2 = nu1, se = se1, se2 = se2) else out <- list(par = theta0, par2 = nu1)
  return(out)
} 
