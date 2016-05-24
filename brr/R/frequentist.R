#' @name FrequentistInference
#' @rdname frequentist 
#' 
#' @title Frequentist inference about the relative risk
#' @description Frequentist confidence intervals about the relative risk: 
#' binomial interval (\code{rr_interval_binomial}) and Sahai and Khurshid confidence interval (\code{rr_interval_SK})
#' 
#' @details The binomial interval (\code{rr_interval_binomial}) is the classical 
#' confidence interval obtained by conditionning on the sum \code{x+y} of 
#' the two counts. The same interval 
#' is implemented in the \code{rateratio.test} package. 
#' The Sahai and Khurshid interval (\code{rr_interval_SK}) is an unconditional confidence interval. See the 
#' reference for more details and a study of its performance.
#' 
#' @references S. Laurent, C. Legrand: A Bayesian framework for the 
#' ratio of two Poisson rates in the context of vaccine efficacy trials. 
#' ESAIM, Probability & Statistics 16 (2012), 375--398.
#' 
#' @param S,T sample sizes
#' @param x,y Observed counts
#' @param conf confidence level
#' 
#' @return \code{rr_interval_binomial} and \code{rr_interval_SK} return 
#' the bounds of the confidence interval 
#' in a vector, \code{rr_intervals} returns a list with the two confidence intervals
#' 
#' @examples 
#' x <- 3; y <- 10; S <- 100; T <- 100
#' rr_intervals(x, y, S, T)
#' brr_intervals(x, y, S, T)
#' 
#' @importFrom stats qnorm qbeta
NULL 
#' 
#' @rdname frequentist
#' @export
rr_interval_SK <- function(x, y, S, T, conf=0.95){
  alpha <- (1-conf)/2
  z <- qnorm(alpha, 0, 1, lower.tail=FALSE)
  D <- x+y+1-0.25*z^2
  if(!(x==0 & y==0) & D>0){
    t1 <- sqrt((x+0.5)*(y+0.5))
    t2 <- 0.5*z*sqrt(D)
    d <- y+0.5-0.25*z^2
    LB <- ((t1-t2)/d)^2
    UB <- ((t1+t2)/d)^2		
  }else{
    LB<-0 
    UB <- Inf
  }
  bounds <- T/S*c(LB,UB)
  return(bounds)
}
#' 
#' @rdname frequentist
#' @export
rr_interval_binomial <- function(x, y, S, T, conf=0.95){
  alpha <- (1-conf)/2
  n <- x+y
  p.L <- function(alpha) {
    if (x == 0) 
      0
    else qbeta(alpha, x, n - x + 1)
  }
  p.U <- function(alpha) {
    if (x == n) 
      1
    else qbeta(1 - alpha, x + 1, n - x)
  }
  L <- p.L(alpha)
  U <- p.U(alpha)
  bounds <- T/S*c(L,U)/(1-c(L,U))
  return(bounds)
}
#' 
#' @rdname frequentist
#' @export
rr_intervals <- function(x, y, S, T, conf=0.95){
  out <- list(binomial=NULL, SK=NULL)
  out$binomial <- rr_interval_binomial(x, y, S, T, conf)
  out$SK <- rr_interval_SK(x, y, S, T, conf)
  attr(out, "level") <- conf
  return(out)
}
