# approximation to derivative of MVN rectangle probabilities using methods
# in Joe (1995), JASA

# lb=lower limit, ub=upper limit, 
# mu = mean vector
# sigma= covariance matrix (this is assumed to be positive definite)
#    matrix functions in R/Splus can be used to check positive definiteness
# k=argument of lb/ub for deriv, ksign=1 for upper, -1 for lower
# eps is tolerance for integration in mvn approx
# nsim is used for m>=7 if random permutations are used in the
#   approximation method
mvn.deriv.margin<- function(lb, ub, mu, sigma, k, ksign, type=1, eps = 1.e-05, nsim=0)
{
  #if(type==1)
  #{ if(!is.loaded(symbol.C("mvndu"))) dyn.load("./mvn.so") }
  #else
  #{ if(!is.loaded(symbol.C("mvndu2"))) dyn.load("./mvn.so") }
  m <- length(ub)
  if(m>=8 && nsim==0) nsim<-10000
  if(m!=length(lb)) stop("lengths of w and x must be the same")
  if(k<1 || k>m) stop("k must be between 1 and dim(lb)")
  if(abs(ksign)!=1) stop("ksign is -1 or 1 for lower/upper limit")
  tem <- sqrt(diag(sigma))
  w <- (lb - mu)/tem
  x <- (ub - mu)/tem
  tem <- diag(1/tem)
  corr <- tem %*% sigma %*% tem
  corr<-c(corr)
  if(eps<1.e-06) eps<- 1.e-06
  
  #print(w)
  #print(x)
  #print(corr)
  ks=k*ksign
  #print(ks)
  if(type==1)
  { out <- .C("r_mvndu",
      as.integer(m), as.double(w), as.double(x), as.double(corr), 
      as.integer(ks), as.integer(nsim),
      as.double(eps), ifail = as.integer(0), deriv=as.double(0))
  }
  else
  { out <- .C("r_mvndu2",
      as.integer(m), as.double(w), as.double(x), as.double(corr), 
      as.integer(ks), as.integer(nsim),
      as.double(eps), ifail = as.integer(0), deriv=as.double(0))
  }
  # need to modify derivative when sigma is not correlation
  #print(out$deriv)
  der=out$deriv
  der=der/sqrt(sigma[k,k])
  out <- list(deriv = der, ifail = out$ifail)
  out
}


# lb=lower limit, ub=upper limit
# mu = mean vector
# sigma= covariance matrix (this is assumed to be positive definite)
#    matrix functions in R/Splus can be used to check positive definiteness
# (j1, k1) index of correlation matrix for derivative
# this is not derivative wrt a covariance!!
# eps is tolerance for integration in mvn approx
# nsim is used for m>=7 if random permutations are used in the
#   approximation method
mvn.deriv.rho<- function(lb, ub, mu, sigma, j1, k1, type=1, eps = 1.e-05, nsim=0)
{
  #if(type==1)
  #{ if(!is.loaded(symbol.C("mvndrh"))) dyn.load("./mvn.so") }
  #else
  #{ if(!is.loaded(symbol.C("mvndrh2"))) dyn.load("./mvn.so") }
  m <- length(ub)
  if(m>=8 && nsim==0) nsim<-10000
  if(m!=length(lb)) stop("lengths of w and x must be the same")
  if(j1<1 || j1>m) stop("k1 must be between 1 and dim(lb)")
  if(k1<1 || k1>m) stop("k1 must be between 1 and dim(lb)")
  if(j1==k1) stop("j1 and k1 should be different indices")
  tem <- sqrt(diag(sigma))
  w <- (lb - mu)/tem
  x <- (ub - mu)/tem
  tem <- diag(1/tem)
  corr <- tem %*% sigma %*% tem
  corr<-c(corr)
  if(eps<1.e-06) eps<- 1.e-06
  
  #print(w)
  #print(x)
  #print(corr)
  if(type==1)
  { out <- .C("r_mvndrh",
      as.integer(m), as.double(w), as.double(x), as.double(corr), 
      as.integer(j1), as.integer(k1), as.integer(nsim),
      as.double(eps), ifail = as.integer(0), deriv=as.double(0))
  }
  else
  { out <- .C("r_mvndrh2",
      as.integer(m), as.double(w), as.double(x), as.double(corr), 
      as.integer(j1), as.integer(k1), as.integer(nsim),
      as.double(eps), ifail = as.integer(0), deriv=as.double(0))
  }
  # need to modify derivative when sigma is not correlation
  out <- list(deriv = out$deriv, ifail = out$ifail)
  out
}

