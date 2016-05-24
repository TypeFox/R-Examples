#' @name PPS
#' @aliases dPPS pPPS hPPS qPPS rPPS
#' @title The Pareto Positive Stable (PPS) distribution
#' @description Density, distribution function, hazard function, quantile function and random generation for the Pareto Positive Stable (PPS) distribution with parameters \code{lam}, \code{sc} and \code{v}.
#' @param x vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of random values to return.
#' @param lam vector of (non-negative) first shape parameters.
#' @param sc vector of (non-negative) scale parameters.
#' @param v vector of (non-negative) second shape parameters.
#' @param log logical; if TRUE, probabilities/densities p are returned as \eqn{log(p)}.
#' @param log.p logical; if TRUE, probabilities/densities p are returned as \eqn{log(p)}.
#' @param lower.tail logical; if TRUE (default), probabilities are \eqn{P[X \le x]}, otherwise, \eqn{P[X > x]}.
#' @details 
#' The PPS distribution has density
#' \deqn{f(x) = \lambda \nu [log(x / \sigma)] ^ (\nu-1) exp(- \lambda [log(x / \sigma)] ^ \nu) / x,}
#' cumulative distribution function
#' \deqn{F(x) = 1 - exp(- \lambda [log(x / \sigma) ^ \nu]),}
#' quantile function
#' \deqn{Q(p) = \sigma exp([- (1 / \lambda) log(1 - p)] ^ (1 / \nu))}
#' and hazard function
#' \deqn{\lambda \nu (log(x / \sigma)) ^ (\nu - 1)  x ^ (-1).}
#' See Sarabia and Prieto (2009) for the details about the numbers random generation.
#' @return
#' \code{dPPS} gives the (log) density, \code{pPPS} gives the (log) distribution function, \code{qPPS} gives the quantile function, and \code{rpois} generates random samples.   
#' Invalid parameters will result in return value \code{NaN}, with a warning. 
#' The length of the result is determined by \code{n} for \code{rPPS}, and is the common length of the numerical arguments for the other functions.
#' @references Sarabia, J.M and Prieto, F. (2009). The Pareto-positive stable distribution: A new descriptive model for city size data, \emph{Physica A: Statistical Mechanics and its Applications}, \bold{388}(19), 4179-4191.
#' @examples
#' print(x <- sort(rPPS(10, 1.2, 100, 2.3)))
#' dPPS(x, 1.2, 100, 2.3)
#' pPPS(x, 1.2, 100, 2.3)
#' qPPS(pPPS(x, 1.2, 100, 2.3), 1.2, 100, 2.3)
#' hPPS(x, 1.2, 100, 2.3)

#' @rdname dPPS
#' @export
dPPS <-
  function(x,lam,sc,v, log = FALSE){
    if (lam<=0 | sc<=0 | v<=0) salida<-0
    else{
      salida<-rep(NA,length(x))
      salida[(x<sc)==1]<-rep(0,sum(x<sc))
      salida[(x<sc)==0]<-lam*v*(log(x[(x<sc)==0]/sc))^(v-1)*x[(x<sc)==0]^(-1)*exp(-lam*(log(x[(x<sc)==0]/sc))^v)
    }
    if (log == TRUE) salida <- log(salida)
    return(salida)
  }

#' @rdname dPPS
#' @export
hPPS <-
  function(x,lam,sc,v){
    if (lam<=0 | sc<=0 | v<=0) salida<-0
    else{
      salida<-rep(NA,length(x))
      salida[(x<sc)==1]<-rep(0,sum(x<sc))
      salida[(x<sc)==0]<-lam*v*(log(x[(x<sc)==0]/sc))^(v-1)*x[(x<sc)==0]^(-1)
    }
    return(salida)
  }

#' @rdname dPPS
#' @export
pPPS <-
  function(x,lam,sc,v, lower.tail = TRUE, log.p = FALSE){
    if (lam<=0 | sc<=0 | v<=0) salida<-0
    else{
      salida<-rep(NA,length(x))
      salida[(x<sc)==1]<-rep(0,sum(x<sc))
      salida[(x<sc)==0]<-1-exp(-lam*(log(x[(x<sc)==0]/sc))^v)
    }
    if (lower.tail == FALSE) salida <- 1 - salida
    if (log.p == TRUE) salida <- log(salida)
    return(salida)
  }

#' @rdname dPPS
#' @export
qPPS <-
  function(p, lam, sc, v, lower.tail = TRUE, log.p = FALSE){
    if (lower.tail == FALSE) p <- 1 - p
    if (log.p == TRUE) p <- exp(p)
    if (lam <= 0 | sc <= 0 | v <= 0) salida <- 0
    else{
      salida <- rep(NA, length(p))
      salida[((p<0)|(p>1))==1] <- rep(0, sum(((p<0)|(p>1))))
      salida[((p<0)|(p>1))==0] <- sc * exp((-(1/lam) * log(1-p[((p<0)|(p>1))==0])) ^ (1/v))
    }
    return(salida)
  }


#' @importFrom stats rweibull
#' @rdname dPPS
#' @export
rPPS <-
  function(n,lam,sc,v){
    weib <- rweibull(n, shape = v, scale = 1)
    sc*exp(lam^(-1/v)*weib)
  }