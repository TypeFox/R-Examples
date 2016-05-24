

greek_utf8 <- function(letter){ #' Unicode encoding of Greek letters
  if(letter %in% c("mu", "phi", "lambda")){
    return(switch(letter, 
                  mu="\u03BC", #rawToChar(as.raw(c(0xc2, 0xb5))) (micron)
                  phi= "\u03d5", # rawToChar(as.raw(c(0xcf, 0x95))), # intToUtf8(0x03D5L)
                  lambda="\u03BB"))
  }else{
    return(letter)
  }
}


#' @importFrom gsl hyperg_2F1
Gauss2F1 <- function(a, b, c, x){ #' Gauss hypergeometric function
  if(x>=0 & x<1){
    hyperg_2F1(a,b,c,x)
  }else{
    #exp(log(hyperg_2F1(c-a,b,c,1-1/(1-x))) - b*log1p(-x))
    hyperg_2F1(c-a,b,c,1-1/(1-x)) / (1-x)^b
  }
}

# Inverse cdf of a discrete distribution
# 
# @param pmf a probability mass function
# @param p probability
# @param ... arguments passed to \code{pmf}
# 
# @examples
# icdf(dpois, 0.5, lambda=10)
# qpois(0.5, 10)
# 
# @export 
icdf <- function(pmf, p, ...){
  q <- 0
  prob <- pmf(0, ...)
  while(prob < p){
    q <- q+1
    prob <- prob + pmf(q, ...)
  }
  return(q)
}

# Moment of a discrete distribution
# 
# @param pmf a probability mass function
# @param k order
# @param accuracy accuracy
# @param ... arguments passed to \code{pmf}
# 
# @examples
# dd_moment(dpois, lambda=5)
# dd_moment(dpois, lambda=5.5795791557050280)
# dd_moment(dpois, k=2, lambda=5)
# @export
dd_moment <- function(pmf, k=1, accuracy=.Machine$double.eps, ...){
 m0 <- 0
 x <- 1
 px <- pmf(x, ...)
 m1 <- m0+x^k*px
 while(px==0 || (m1-m0)>accuracy){
   x <- x + 1
   px <- pmf(x, ...)
   m0 <- m1
   m1 <- m0+x^k*px
 }
 return(m1)
}

# Integration range for beta expectation of a nonnegative function on (0,1)
# @param c,d shape parameters
# @param f the function 
#' @importFrom stats dbeta qbeta
beta_integration_range <- function(c, d, f, min=100, accuracy=1e-8){
  # 
  if(c+d<min || (c<=1 && d<=1)) return(c(0,1))
  if(c>1 && d>1) accuracy <- accuracy/2
  if(d>1){
    i <- -3
    ubound <- qbeta(1-10^i, c, d)
    if(ubound>.999){
      ubound <- 1
    }else{
      value <- f(ubound)*dbeta(ubound,c,d)
      while(value>accuracy){
        i <- i-1
        ubound <- qbeta(1-10^i, c, d)
        value <- f(ubound)*dbeta(ubound,c,d)
      }
    }
  }
  if(c<=1) return(c(0,ubound)) 
  if(d<=1) return(rev(1-beta_integration_range(d, c, function(x) f(1-x), accuracy)))
  #
  i <- -3
  lbound <- qbeta(10^i, c, d)
  if(lbound<0.001){
    lbound <- 0
  }else{
    value <- f(lbound)*dbeta(lbound,c,d)
    while(value>accuracy){
      i <- i-1
      lbound <- qbeta(10^i, c, d)
      value <- f(lbound)*dbeta(lbound,c,d)
    }
  }
  return(c(lbound,ubound))
}
