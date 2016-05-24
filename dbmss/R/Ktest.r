Ktest <- function(X, r)
{

  # Verify X
  if (!inherits(X, "ppp"))
    stop("X is not of class ppp") 
  # Verify r
  if (!is.numeric(r) && !is.vector(r)) 
    stop("r must be a numeric vector")
  if (length(r) < 2) 
    stop(paste("r has length", length(r), "- must be at least 2"))
  if (any(diff(r) <= 0)) 
    stop("successive values of r must be increasing")  
  
  # Utilities
  # r : radius
  # w, l : width and length of the window
  #-------------------------------------------------------------------------------
  # ern : Probability for a point to be in another's neighborhood.(pi * r^2 / (w*l) corrected from edge effects)
  ern <- function(r, w, l)
  {
    return(r*r/w/l*(pi-4*r/3*(1/w+1/l)+r*r/2/w/l))
  }
  #-------------------------------------------------------------------------------
  # espKc : Expectation of K, when density is known
  espKc <- function(r,w,l)
  {
    return(w*l*ern(r,w,l))
  }
  #-------------------------------------------------------------------------------
  # espKi : Expectation of K, when density is unknown. rho*w*l is estimated by the number of points (N).
  #        Difference between espKi and espKc goes to 0 when points are more than 20.
  espKi <- function(r, N, w, l)
  {
    return(espKc(r, w, l)*(1-(1+N)*exp(-N)))
  }
  #-------------------------------------------------------------------------------
  # eh : Expectation of h_1^2(U, r)
  eh <- function(r, w, l)
  {
    return(r^5*(w+l)/2/w^3/l^3*(8*pi/3-256/45) +r^6/w^3/l^3*(11*pi/48+8/9-16/9*(w+l)*(w+l)/w/l) +4/3*r^7*(w+l)/w^4/l^4 -r^8/4/w^4/l^4)
  }
  # ------------------------------------------------------------------------------
  # sigmaKi : Variance matrix for the estimator of K (unknown rho), normalized by multiplication by par w*l*rho
  # vec     : Vector of distances to compute K (example : c(1, 2) to calculate K(1) and K(2)
  # N       : number of points
  # w       : width of the rectangle window
  # l       : length of the rectangle window
  sigmaKi <- function(vec, N, w, l)
  {
    # Intermediate computations
    d <- length(vec)
    ern_ <- ern(vec, w, l)
    
    c1 <- 2*w*w*l*l/(N*(N-1))
    c2 <- w*w*l*l*exp(-N)*(1+N)*(1-exp(-N)-N*exp(-N))
    c3 <- 2*(N-2)*c1
  
    # Preparation of a square matrix
    sigmaKi_ <- matrix(nrow = d, ncol = d)
  
    # Top half of the matrix
    for (i in 1:(d-1))
    {
      for (j in i:d)
      {
        covh1_ <- covh1(vec[i], vec[j], w, l)
        sigmaKi_[i, j] <- c1*ern_[min(i, j)]+(c2-c1)*ern_[i]*ern_[j]+c3*covh1_
      }
    }
    # Diagonal
    for (i in 1:d)
    {
      sigmaKi_[i, i] <- c1*ern_[i]+(c2-c1)*ern_[i]*ern_[i]+c3*eh(vec[i], w, l)
    }
    # Bottom half
    for (j in 1:(d-1))
    {
      for (i in j:d)
      {
        sigmaKi_[i, j] <- sigmaKi_[j, i]
      }
    }
    # Result                                                            
    return(sigmaKi_)
  }
  #-------------------------------------------------------------------------------
  # invsqrtmat : Transforms a matrix into the square root of its inverse
  #             such as invsqrtmat %*% mat %*% t(invsqrtmat) = Id
  invsqrtmat <- function(mat)
  {
    if(length(mat)>1)
    {
      e <- eigen(mat)
      # Eigen vectors
      p <- e$vectors
      # Square roots of eigen values
      d <- sqrt(e$values)
      # put into a diagonal matrix
      rd <- diag(d)
      # Resolution
      return(solve(p%*%rd))
    }
    else
    {
      return(1/sqrt(mat))
    }
  }
  #-------------------------------------------------------------------------------
  # covh1  : covariance of h1
  # r1, r2 : distances
  # w      : width of the rectangle window
  # l      : length of the rectangle window
  covh1 <- function(r1, r2, w, l)
  {
    # Values of the product in the corner. Coordinates are x'=(n-xi)/r2 without normalization in r'^2r^2/w^2/l^2
    # This function must be inside covh1 because it needs to use local variables, that can not be passed as parameters because adaptIntegrate forbides it.
    corner <- function(x){(foncA1(ra2*x/ra1, ra1, w, l)+foncA2(ra2*x/ra1, ra1, w, l)+foncA3(ra2*x/ra1, ra1, w, l)+foncA4(ra2*x/ra1, ra1, w, l))*(foncA3(x, ra2, w, l)+foncA4(x, ra2, w, l))}
    
    # r values must be ordered
    ra1 <- min(r1, r2)
    ra2 <- max(r1, r2)
  
    # Size ratios parameters
    rapr <- ra1/ra2
    r12 <- ra1*ra1/w/l
    r22 <- ra2*ra2/w/l
  
    # Biases, normalized by n^2/r^2
    b1 <- brn(ra1, w, l)
    b2 <- brn(ra2, w, l)
    
    # Numerical computing of the elliptic integral
    int2 <- stats::integrate(integrand3, lower=0, upper=1, r1=ra1, r2=ra2)
    intcorner <- cubature::adaptIntegrate(corner, lowerLimit=c(0, 0), upperLimit=c(1, 1))
  
    # line 1
    covh1_ <- (w-2*r2)/w*(l-2*r2)/l*b1*b2
    # line 2
    covh1_ <- covh1_+2*(w+l-4*r2)*r2/w/l*b1*(b2-foncG(1))
    # line 3
    covh1_ <- covh1_+2*(w+l-4*r2)*r1/w/l*(int2$value-b2*foncG(1))+4*r22*intcorner$integral
   
    # multiplication by the common factor
    covh1_ <- covh1_*r12*r22
  
    # Result                                                            
    return(covh1_)
  }
  
  foncg   <- function(x){if (x<=1) { return(acos(x)-x*sqrt(1-x*x))} else {return(0)}}
  foncg01 <- function(x){acos(x)-x*sqrt(1-x*x)} 
  foncG   <- function(x){x*acos(x)-sqrt(1-x*x)*(2+x*x)/3+2/3}
  
  # Values of h1 on different zones without normalization by r^2/w/l
  brn    <- function(r, w, l){4*r*(w+l)/(3*w*l)-r*r/(2*w*l)}
  indic  <- function(a){as.numeric(a)}
  foncA1 <- function(x, r, w, l){brn(r, w, l)*indic(x[1] >= 1)*indic(x[2] >= 1)}
  foncA2 <- function(x, r, w, l){(brn(r, w, l)-foncg(x[2]))*indic(x[1] >= 1)*indic(x[2] < 1)+(brn(r, w, l)-foncg(x[1]))*indic(x[2] >= 1)*indic(x[1] < 1)}
  foncA3 <- function(x, r, w, l){(brn(r, w, l)-foncg(x[1])-foncg(x[2]))*indic(x[1] < 1)*indic(x[2] < 1)*indic(x[1]^2+x[2]^2 > 1)}
  foncA4 <- function(x, r, w, l){(brn(r, w, l)+x[1]*x[2]-(foncg(x[1])+foncg(x[2]))/2-pi/4)*indic(x[1] < 1)*indic(x[2] < 1)*indic(x[1]^2+x[2]^2 <= 1)}
  
  integrand3 <- function(x, r1, r2){foncg01(r1*x/r2)*foncg01(x)}
  #-------------------------------------------------------------------------------
  # End of utilities
  
  # fails if the window is not a rectangle
  if (!is.rectangle(X$window)) stop("The window must be a rectangle to apply Ktest.")
    
  #  value = NA if the number of points is less than 2  
  if (X$n > 1)
    {
    # Width - Length
    l <- X$window$xrange[2]-X$window$xrange[1]
    w <- X$window$yrange[2]-X$window$yrange[1]
    # espKi : expectation of K, rho unknown, calculated according to the number of points.
    espKi_ <- espKi(r, X$n, w, l)

    # Computed variance. Depends on the the observed number of points and the window.
    sigmaKi_ <- sigmaKi(r, X$n, w, l)
   
    # Distance matrix.
    pairdist_ <- pairdist.ppp(X)  # Requires a lot of RAM. Limited to 8, 000 points with 32-bit versions of R.

    # Estimation of K
    # NbPairs : number of pairs of points less than r apart (it's a vector, one value for each r)
    # pairdist_ >0 eliminates distance from a point to itself
    # *1.0 is to convert values to real and avoid infinite values in NbPairs*w*l later
    NbPairs <- sapply(r, function(d) sum(pairdist_ > 0 & pairdist_ < d)*1.0)
    # Kest gets the estimator of K, centered on the expected value.
    Kest <- NbPairs*w*l/(X$n*(X$n-1))-espKi_
    
    TestVector <- invsqrtmat(sigmaKi_) %*% Kest 
    return(1 - stats::pchisq(sum(TestVector*TestVector), length(r)))
  }
  else return(NA) 
}
