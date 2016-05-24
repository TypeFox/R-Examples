findRoots <- function(coefs,k=1) {
  # Return complex roots of polynomial equation
  #   coefs[,3]*x^(2*k) + coefs[,2]*x^k + coefs[,1] = 0

  # Recode coefs as a matrix
  attr(coefs,"dim") <- c(length(coefs)/3,3)

  # Sanity check
  if (any(coefs[,3]==0)) stop("Leading coefficients must be non-zero")
  if (any((coefs[,1]==0) & (coefs[,2]==0))) stop("Double roots not acceptable in present implementation. I suggest that you add a contant term to the L-operator.")
  
  # Use method described in Press et al, "Numerical Recipies in C",
  # second edition, Section 5.6.
  q <- -(coefs[,2]+(1-2*(coefs[,2]<0))*
         sqrt(complex(real=coefs[,2]^2-4*coefs[,3]*coefs[,1])))/2
  x <- cbind(q/coefs[,3],coefs[,1]/q)

  # Build all roots
  x <- kronecker(x^(1/k),matrix(complex(argument=2*pi*(1:k)/k),1,k))
  x <- t(apply(x,1,sort))

  # Return result
  return(list(left=as.matrix(x[,k:1]),right=as.matrix(x[,k+1:k])))
}
