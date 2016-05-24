#' Analytically integrate logistic detection function
#'
#' Computes integral (analytically) over x from 0 to width of a logistic
#' detection function; For reference see integral #526 in CRC Std Math
#' Table 24th ed
#'
#' @param x matrix of data
#' @param models list of model formulae
#' @param beta parameters of logistic detection function
#' @param width transect half-width
#'
#' @author Jeff Laake
#'
integratelogistic.analytic <- function(x, models, beta, width){

  # integral of f(x) = 1/(a+b*exp(px)) - in this case a=1 and b=exp(a) where:
  #  * a is the part of the logistic that doesn't involve distance (x).
  #  * b is the portion of the integral that doesn't include distance which
  #    is computed by setting x=0

  b <- logit(logisticbyz(x,0,models,beta))

  # The coefficient for distance is obtained by computing
  # logit(f(x=0))-logit(f(x=1)) = a + p*0 - (a+p*1) = = -p.  The sign is changed
  # so the logistic 1/(1+b*exp(-px)) matches the form of the function
  p <- b - logit(logisticbyz(x,1,models,beta))

  # what is needed is exp(-b) for the logistic form
  b <- exp(-b)

  # if p=0 then the function is constant with respect to x (distance) so the
  # integral is straightfoward
  int <- width/(1+b)

  # For those values with p !=0, the following is the integral for 0 to width
  int[p!=0] <- width - log(1+b[p!=0]*exp(p[p!=0]*width))/p[p!=0] +
                 log(1+b[p!=0])/p[p!=0]

  return(int)
}
