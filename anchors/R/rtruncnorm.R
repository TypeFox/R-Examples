rtruncnorm <- function(n,lower,upper) {

  if (n > 1 && length(lower)==1) {
    lower <- rep(lower,n)
  } else {
    stop("Length of lower bound neither scalar nor equal to n\n")
  }
  if (n > 1 && length(upper)==1) {
    upper <- rep(upper,n)
  } else {
    stop("Length of upper bound neither scalar nor equal to n\n")
  }

  .C("rtruncnorm",
                n = as.integer(n),
                xlower=as.double(lower),
                xupper=as.double(upper),
                rval=as.double(rep(-99999,n)), PACKAGE="anchors")$rval
}

rtruncnorm1 <- function(n,mu,sigma,lower,upper) {

  if (n >= 1 && n%%length(lower)==0 && n%%length(upper)==0 && n%%length(mu)==0  && n%%length(sigma)==0 ) {
    lower <- rep(lower,n/length(lower))
    upper <- rep(upper,n/length(upper))
    mu    <- rep(mu   ,n/length(mu))
    sigma <- rep(sigma,n/length(sigma))
  } else {
    stop("Lengths not equal\n")
  }  

  .C("rtruncnorm1",
     n = as.integer(n),
     mu = as.double(mu),
     sigma = as.double(sigma),
     xlower=as.double(lower),
     xupper=as.double(upper),
     PACKAGE="anchors")$mu
}

rtruncnorm2 <- function(n,mu,sigma,lower,upper) {

  if (n >= 1 && n%%length(lower)==0 && n%%length(upper)==0 && n%%length(mu)==0  && n%%length(sigma)==0 ) {
    lower <- rep(lower,n/length(lower))
    upper <- rep(upper,n/length(upper))
    mu    <- rep(mu   ,n/length(mu))
    sigma <- rep(sigma,n/length(sigma))
  } else {
    stop("Lengths not equal\n")
  }  

  .C("rtruncnorm2",
     n = as.integer(n),
     mu = as.double(mu),
     sigma = as.double(sigma),
     xlower=as.double(lower),
     xupper=as.double(upper),
     PACKAGE="anchors")$mu
}

convolve2 <- function(a, b)
  .C("convolve2",
     as.double(a),
     as.integer(length(a)),
     as.double(b),
     as.integer(length(b)),
     ab = double(length(a) + length(b) - 1),PACKAGE="anchors")$ab

runif9 <- function(a,b)
  .C("runif9",
     as.double(a),
     as.double(b),
     ab = as.double(-99),
     PACKAGE="anchors")$ab

