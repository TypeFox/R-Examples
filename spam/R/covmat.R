# This is file ../spam/R/covmat.R
# This file is part of the spam package, 
#      http://www.math.uzh.ch/furrer/software/spam/
# by Reinhard Furrer [aut, cre], Florian Gerber [ctb]
     


# construct various precision matrices



covmat <- function(h, theta, ... , type="sph") {
  avtype <- c("exponential", "spherical", "nugget",
              "wu1","wu2","wu3","wendland1","wendland2", "matern")
  method <- pmatch(tolower(type), avtype)
  if (is.na(method)) 
     stop("Covariance function not implemented yet. Please ask for.")
  switch(method,
         return(cov.exp(h, theta, ...)),
         return(cov.sph(h, theta, ...)),
         return(cov.nug(h, theta, ...)),
         return(cov.wu1(h, theta, ...)),
         return(cov.wu2(h, theta, ...)),
         return(cov.wu3(h, theta, ...)),
         return(cov.wend1(h, theta, ...)),
         return(cov.wend2(h, theta, ...)),
         return(cov.mat(h, theta, ...)))
}

.par.check.cov <- function(theta,nr=2){
  if (any(theta<0)) {
    warning('Parameters coerced to positive values')
    theta <- abs(theta)
  }
  nt <- length(theta)
  
  if (nt < nr) 
    return( c( theta, rep(1, nr-nt), 0))
  return( c( theta, 0)[1:(nr+1)])
}
  

cov.sph <- function(h, theta, ... , eps= .Spam$eps) {
  theta <- .par.check.cov(theta)
  
  if (is.spam(h)) {
    tmp <- h@entries/theta[1]
    h@entries <- ifelse(tmp < eps, theta[2] + theta[3],
                        ifelse(tmp < 1, theta[2] * (1 - 1.5 * tmp + 0.5 * tmp^3), 0))
    return( h)    
  } else {
    h <- h/theta[1]
    ifelse(h < eps, theta[2] + theta[3],
           ifelse(h < 1, theta[2] * (1 - 1.5 * h + 0.5 * h^3), 0))
  }
  
    
}
cov.wend1 <- function(h, theta,  ... , eps= .Spam$eps) {
  # is \phi_{3,1} in the 98 paper and \psi_{3,1} in the 95 paper
  theta <- .par.check.cov(theta)
  
  if (is.spam(h)) {
    tmp <- h@entries/theta[1]
    h@entries <- ifelse(tmp < eps, theta[2] + theta[3],
                        ifelse(tmp < 1,
                               theta[2]  * ((1 - tmp)^4*(4*tmp+1)), 0))
    return( h)    
  } else {
    h <- h/theta[1]
    ifelse(h < eps, theta[2] + theta[3],
                        ifelse(h < 1,
                               theta[2]  * ((1 - h)^4*(4*h+1)), 0))

  }
}

cov.wend2 <- function(h, theta,  ... , eps= .Spam$eps) {
  # is \phi_{3,2} in the 98 paper and \psi_{4,2} in the 95 paper
   theta <- .par.check.cov(theta)
  
  if (is.spam(h)) {
    tmp <- h@entries/theta[1]
    h@entries <- ifelse(tmp < eps, theta[2] + theta[3],
                        ifelse(tmp < 1,
                               theta[2] * ((1 - tmp)^6*(35*tmp^2+18*tmp+3))/3, 0))
    return( h)    
  } else {
    h <- h/theta[1]
    ifelse(h < eps, theta[2] + theta[3],
                        ifelse(h < 1,
                               theta[2] * ((1 - h)^6*(35*h^2+18*h+3))/3, 0))

  }
}
cov.wu1 <- function(h, theta, ... ,  eps= .Spam$eps) {
   theta <- .par.check.cov(theta)
  
  if (is.spam(h)) {
    tmp <- h@entries/theta[1]
    h@entries <- ifelse(tmp < eps, theta[2] + theta[3],
                        ifelse(tmp < 1,
                               theta[2] * ((1 - tmp)^3*(1+3*tmp+tmp^2)), 0))
    return( h)    
  } else {
    h <- h/theta[1]
    ifelse(h < eps, theta[2] + theta[3],
                        ifelse(h < 1,
                               theta[2] * ((1 - h)^3*(1+3*h+h^2)), 0))

  }
}
cov.wu2 <- function(h, theta,  ... , eps= .Spam$eps) {
   theta <- .par.check.cov(theta)
  
  if (is.spam(h)) {
    tmp <- h@entries/theta[1]
    h@entries <- ifelse(tmp < eps, theta[2] + theta[3],
                        ifelse(tmp < 1,
                               theta[2] * ((1 - tmp)^4*(4+16*tmp+12*tmp^2+3*tmp^3))/4, 0))
    return( h)    
  } else {
    h <- h/theta[1]
    ifelse(h < eps, theta[2] + theta[3],
                        ifelse(h < 1,
                               theta[2] * ((1 - h)^4*(4+16*h+12*h^2+3*h^3))/4, 0))

  }
}
cov.wu3 <- function(h, theta,  ... , eps= .Spam$eps) {
   theta <- .par.check.cov(theta)
  
  if (is.spam(h)) {
    tmp <- h@entries/theta[1]
    h@entries <- ifelse(tmp < eps, theta[2] + theta[3],
                        ifelse(tmp < 1,
                               theta[2] * ((1 - tmp)^6*(1+6*tmp+41/3*tmp^2+12*tmp^3+5*tmp^4+5/6*tmp^5)), 0))
    return( h)    
  } else {
    h <- h/theta[1]
    ifelse(h < eps, theta[2] + theta[3],
                        ifelse(h < 1,
                               theta[2] * ((1 - h)^6*(1+6*h+41/3*h^2+12*h^3+5*h^4+5/6*h^5)), 0))

  }
}



cov.mat <- function(h, theta,  ... ,  eps= .Spam$eps)
{
  theta <- .par.check.cov(theta,3)
  if (is.spam(h)) {
    tmp <- h@entries/theta[1]
    h@entries <- ifelse(tmp < eps, theta[2] + theta[4],
                        theta[2] * (((2^(-(theta[3] - 1)))/gamma(theta[3])) * 
                                    (tmp^theta[3]) * besselK(tmp, nu = theta[3])))
       
    return( h)    
  } else {
    h <- h/theta[1]
    ifelse(h < eps, theta[2] + theta[4],
                        theta[2] * (((2^(-(theta[3] - 1)))/gamma(theta[3])) * 
                                    (h^theta[3]) * besselK(h, nu = theta[3])))
  }
}



cov.exp <- function(h, theta,  ... ,  eps= .Spam$eps)
{
  theta <- .par.check.cov(theta,2)
  if (is.spam(h)) {
    tmp <- h@entries/theta[1]
    h@entries <- ifelse(tmp < eps, theta[2] + theta[3],
                        theta[2] * exp( -tmp))
    return( h)    
  } else {
    h <- h/theta[1]
    ifelse(h < eps, theta[2] + theta[3],
           theta[2] * exp( -h))
  }
}

cov.nug <- function(h, theta,  ... , eps= .Spam$eps)
{
  theta <- .par.check.cov(theta,0)
  if (is.spam(h)) {
    h@entries <- ifelse(h@entries < eps, theta[1], 0)
    return( h)    
  } else {
    ifelse(h < eps, theta[1], 0)
  }
}

