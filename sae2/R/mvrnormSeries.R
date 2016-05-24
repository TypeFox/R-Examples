mvrnormSeries <- function(NV=1, D, T, sigma.e, 
   rho.dyn, sigma.v.dyn, sigma.u.dyn, rho.u.dyn,
   rho.RY, sigma.v.RY, sigma.u.RY, rho.u.RY,
   tol=1e-6) {
  if (missing(T)) stop("'T' must be specified")
  if (T < 3) stop("'T' must be >2")
  if (class(sigma.e) == "list") {
    if (missing(D)) { 
      D <- length(sigma.e) 
    } else {
      if (D != length(sigma.e)) 
        stop("the length of 'sigma.e' must agree with D")
    }
    lapply(sigma.e, FUN = function(x) { if (!all(dim(x) == c(NV*T, NV*T) )) 
      stop(paste("each element of 'sigma.e' must be a square matrix",
        "with 'NV*T' rows"))
     } )
  } else {
    if (missing(D)) stop("'D' must be specified")
    if (class(sigma.e) != "matrix") 
      stop("'sigma.e' must be a matrix or a list of matrices")
    if (!all(dim(sigma.e) == c(NV*T, NV*T) ))
      stop("'sigma.e' must be a square matrix with 'NV*T' rows")
  }
  sigma <- matrix(0, NV*T, NV*T)
  mu <- rep(0, NV*T)
  if (!(missing(rho.dyn) | missing(sigma.v.dyn) | missing(sigma.v.dyn))) {
    if (missing(rho.dyn) | missing(sigma.v.dyn) | missing(sigma.v.dyn))
      stop("some parameters missing for dynamic model")
    if (NV ==1) {
      if(length(sigma.v.dyn) != 1 |
         length(sigma.u.dyn) != 1)
        stop("wrong length for 'sigma.v.dyn' or 'sigma.u.dyn'")
      sigma <- sigma.v.dyn*Gamma_v_f(rho.dyn, T)[[1]] +
               sigma.u.dyn*Gamma_u_f(rho.dyn, T)[[1]]
    } else {
      if(missing(rho.u.dyn)) stop("'rho.u.dyn' must be specified")
      len.rho <- (NV*(NV-1))/2
      if(length(rho.u.dyn)==1 & NV>2) rho.u.dyn <- rep(rho.u.dyn, len.rho)
      if(length(rho.u.dyn)!= len.rho) 
        stop(paste0("the length of 'rho.u.dyn' must be 1 or (NV*(NV+1))/2"))
      u_corr <- diag(1, NV)
      for (i in 2:NV) for (j in 1:(i-1)){ 
        u_corr[i,j] <- u_corr[j,i] <- rho.u.dyn[((i-2)*(i-1))/2+j]
      }
      if(length(sigma.v.dyn)==1) sigma.v.dyn <- rep(sigma.v.dyn, NV)
      if(length(sigma.u.dyn)==1) sigma.u.dyn <- rep(sigma.u.dyn, NV)
      if(length(sigma.v.dyn)!= NV) 
        stop(paste0("the length of 'sigma.v.dyn' must be 1 or NV"))
      if(length(sigma.u.dyn)!= NV) 
        stop(paste0("the length of 'sigma.u.dyn' must be 1 or NV"))
      sqrt_u <- sqrt(sigma.u.dyn)
      sqrt_v <- sqrt(sigma.v.dyn)
      sigma  <- 
        ((sqrt_u%*%t(sqrt_u)) * u_corr) %x% Gamma_u_f(rho.dyn, T)[[1]] +
        ((sqrt_v%*%t(sqrt_v)) * u_corr) %x% Gamma_v_f(rho.dyn, T)[[1]]
    }
  }
  if (!(missing(rho.RY) | missing(sigma.v.RY) | missing(sigma.v.RY))) {
    if (missing(rho.RY) | missing(sigma.v.RY) | missing(sigma.v.RY))
      stop("some parameters missing for Rao-Yu model")
    Gamma_v <- matrix(1,nrow=T,ncol=T)
    Gamma_b <- matrix(rep(1:T,times=T),nrow=T, ncol=T)
    Gamma_s <- abs(Gamma_b-t(Gamma_b))
    if(rho.RY > 0) {
     Gamma_u <- rho.RY^(Gamma_s)/(1-rho.RY^2)
    } else {
      Gamma_u <- diag(T)
    }
    if (NV ==1) {
      if(length(sigma.v.RY) != 1 |
         length(sigma.u.RY) != 1)
        stop("wrong length for 'sigma.v.RY' or 'sigma.u.RY'")
      sigma <- sigma + sigma.v.RY*Gamma_v +
               sigma.u.RY*Gamma_u
    } else {
      if(missing(rho.u.RY)) stop("'rho.u.RY' must be specified")
      len.rho <- (NV*(NV-1))/2
      if(length(rho.u.RY)==1 & NV>2) rho.u.RY <- rep(rho.u.RY, len.rho)
      if(length(rho.u.RY)!= len.rho) 
        stop(paste0("the length of 'rho.u.RY' must be 1 or (NV*(NV+1))/2"))
      u_corr <- diag(1, NV)
      for (i in 2:NV) for (j in 1:(i-1)){ 
        u_corr[i,j] <- u_corr[j,i] <- rho.u.RY[((i-2)*(i-1))/2+j]
      }
      if(length(sigma.v.RY)==1) sigma.v.RY <- rep(sigma.v.RY, NV)
      if(length(sigma.u.RY)==1) sigma.u.RY <- rep(sigma.u.RY, NV)
      if(length(sigma.v.RY)!= NV) 
        stop(paste0("the length of 'sigma.v.RY' must be 1 or NV"))
      if(length(sigma.u.RY)!= NV) 
        stop(paste0("the length of 'sigma.u.RY' must be 1 or NV"))
      sqrt_u <- sqrt(sigma.u.RY)
      sqrt_v <- sqrt(sigma.v.RY)
      sigma  <- sigma +
        ((sqrt_u%*%t(sqrt_u)) * u_corr) %x% Gamma_u +
        ((sqrt_v%*%t(sqrt_v)) * u_corr) %x% Gamma_v
    }
  }
  if (NV == 1) {
    y <- rep(0, D*T)
    if(class(sigma.e) == "list") {
      for (i in 1:D) {
        y[((i-1)*T + 1):(i*T)] <-
          mvrnorm(mu=mu, Sigma=sigma+sigma.e[[i]], tol=tol)
      }
    } else {
      for (i in 1:D) {
        y[((i-1)*T + 1):(i*T)] <-
          mvrnorm(mu=mu, Sigma=sigma+sigma.e, tol=tol)
      }
    }
  } else {
    y <- matrix(0, nrow=D*T, ncol=NV)
    if(class(sigma.e) == "list") {
      for (i in 1:D) {
        y.temp <- 
            mvrnorm(mu=mu, Sigma=sigma+sigma.e[[i]], tol=tol)
        for (nv in 1:NV) {
           y[((i-1)*T+1):(i*T), nv] <- y.temp[((nv-1)*T+1):(nv*T)]
        }
      }
    } else {
      for (i in 1:D) {
        y.temp <- 
            mvrnorm(mu=mu, Sigma=sigma+sigma.e, tol=tol)
        for (nv in 1:NV) {
           y[((i-1)*T+1):(i*T), nv] <- y.temp[((nv-1)*T+1):(nv*T)]
        }
      }
    }
  }
  return(y)
}

