
bhnoneq <- function(year = NULL, mlen = NULL, ss = NULL, K = NULL, Linf = NULL, 
               Lc = NULL, nbreaks = NULL, styrs = NULL, stZ = NULL, graph = TRUE) {
  if (is.null(mlen)) stop("mean length vector does not exist")
  if (is.null(year)) stop("year vector does not exist")
  if (is.null(ss)) stop("numbers vector does not exist")
  if (!is.numeric(mlen)) stop("vector is not numeric")
  if (is.null(stZ)) stop("Initial Z vector does not exist")
  if (is.null(K)) stop("K not specified")
  if (is.null(Linf)) stop("Linf not specified")
  if (is.null(Lc)) stop("Lc not specified")
  if (length(mlen) != length(year)) stop("vectors have different lengths")
    
  if (is.null(nbreaks)) stop("Number of mortality breaks not specified")
  if (nbreaks > 0 & is.null(styrs)) stop("Starting guesses for years of mortality breaks not specified")
  if (nbreaks > 0 & length(styrs) != nbreaks) stop("Number of starting guesses for years of mortality breaks not equal to number of mortality breaks")
  if (nbreaks > 0 & length(stZ) != nbreaks+1) stop("Number of starting guesses for Z not equal to number of mortality breaks plus one")
    
  dat <- data.frame(year = year, mlen = mlen, ss = ss)
  sigma<-NULL; Lpred<-NULL;results<-NULL
  return_LL <- function(x, dat, nbreaks, K, Linf, Lc) {
  year <- dat$year
  Lbar <- dat$mlen
  ss <- dat$ss
  
  Z <- x[1:(nbreaks+1)]
  ggyr <- x[(nbreaks + 2):length(x)]
  
  nyr <- sum(!is.na(Lbar))
  count <- length(year)
    
  dy <- array(0, dim = c(nbreaks, count))
  
  for(i in 1:nbreaks) {
    for(j in 1:count) {
      dy[i, j] <- ifelse(ggyr[i] >= j, 0, j - ggyr[i])
    }
  }
  
  a <- array(0, dim = c(nbreaks + 1, count))
  s <- array(0, dim = c(nbreaks + 1, count))
  r <- array(0, dim = c(nbreaks + 1, count))
  
  denom <- rep(0, count)
  numersum <- rep(0, count)
  numerator <- rep(0, count)
  
  for(m in 1:count) {
    
    for(i in 1:(nbreaks + 1)) {
      
      a[i, m] <- 1
      r[i, m] <- 1
      
      if (i < (nbreaks + 1)) s[i, m] <- 1 - exp(-(Z[nbreaks + 2 - i] + K) * dy[nbreaks + 1 - i, m])
      if (i == (nbreaks + 1)) s[i, m] <- 1
      
      for (j in 1:(i - 1)) {
        if (i > 1) {
          a[i, m] <- a[i, m] * exp(-Z[nbreaks + 2 - j] * dy[nbreaks + 1 - j, m])
          r[i, m] <- r[i, m] * exp(-(Z[nbreaks + 2 - j] + K) * dy[nbreaks + 1 - j, m])
        }
      }
      
      if (i <= nbreaks) {
        denom[m] <- denom[m] + a[i, m] * 
          ((1 - exp(-Z[nbreaks + 2 - i] * dy[nbreaks + 1 - i, m]))/Z[nbreaks + 2 - i])
      }
      
      if (i == (nbreaks + 1)) {
        denom[m] <- denom[m] + a[i, m]/Z[nbreaks + 2 - i]
      }
      
      numersum[m] <- numersum[m] + r[i, m] * s[i, m]/(Z[nbreaks + 2 - i] + K)
    }
  }
  
  numerator <- Linf * (denom - (1 - Lc/Linf) * numersum)
  Lpred <<- numerator/denom
  
  sigma <<- sqrt(sum(ss * (Lbar - Lpred)^2, na.rm = TRUE)/nyr)
  
  LL <- -nyr * log(sigma) - sum(ss * (Lbar - Lpred)^2, na.rm = TRUE)/(2 * sigma^2)
  return(-1 * LL)
}
return_LL_eq <- function(x, dat, K, Linf, Lc) {
  
  year <- dat$year
  Lbar <- dat$mlen
  ss <- dat$ss
  
  Z <- x
  
  nyr <- sum(!is.na(Lbar))
  count <- length(year)
  
  Lpred <<- rep(Linf * (1 - (Z/(Z+K))*(1 - Lc/Linf)), count)
  sigma <<- sqrt(sum(ss * (Lbar - Lpred)^2, na.rm = TRUE)/nyr)
  
  LL <- -nyr * log(sigma) - sum(ss * (Lbar - Lpred)^2, na.rm = TRUE)/(2 * sigma^2)
  return(-1 * LL)
  
}
  
  if(nbreaks == 0) {
    
    results <- optim(stZ, return_LL_eq, method = "BFGS", dat = dat, K = K, Linf = Linf, Lc = Lc,
                     control = list(maxit = 1e+06, abstol = 1e-07), hessian = TRUE)
  }
  else {
    g_year <- styrs - year[1] + 1
    stpar <- c(stZ, g_year)
    
    results <- optim(stpar, return_LL, method = "BFGS", dat = dat, nbreaks = nbreaks, K = K, Linf = Linf, Lc = Lc, 
                     control = list(maxit = 1e+06, abstol = 1e-07), hessian = TRUE)
  }
    
  npar <- length(results$par) + 1
  AIC <- 2 * results$value + 2 * npar
    
  var <- solve(results$hessian)
  SE <- sqrt(diag(var))
  
  if(nbreaks == 0) {
    tab <- data.frame(Parameter = c("Z", "sigma", "AIC"), 
                      Estimate = c(results$par, sigma, AIC),
                      StdErr = c(SE, NA, NA))
  }
  else {
    vname <- c(paste("Z", seq(1, nbreaks + 1, 1), sep = ""), 
               paste("Y", seq(1, nbreaks, 1), sep = ""), "sigma", "AIC")
    
    est <- c(results$par)
    est[(nbreaks +2):length(est)] <- est[(nbreaks +2):length(est)] + year[1] - 1
    
    tab <- data.frame(Parameter = vname, Estimate = c(est, sigma, AIC), StdErr = c(SE, NA, NA))
  }
  
  if(graph) {
    par(mfrow = c(1, 2))
    plot(year, mlen, xlab = "Year", ylab = "Mean Length", typ = "b")
    lines(year, Lpred, col = "red")
    plot(year, mlen - Lpred, xlab = "Year", ylab = "Residual", typ = "b")
    abline(h = 0, col = "red")
  }
  
  if(results$convergence == 0) status <- paste("Optimizer converged after",results$counts[1],"function calls.") 
  else status <- "Optimizer did not converge."
  
  if(all(eigen(results$hessian)$values > 0)) hess.status <- "Hessian matrix is positive definite."
  else hess.status <- "Hessian matrix is not positive definite. Check parameter estimates."
  
  return(list(summary = tab, convergence = status, hessian = hess.status, 
              results = data.frame(year = year, observed = mlen, predicted = Lpred, residual = mlen - Lpred)))
}
