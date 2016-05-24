
#' @name Post_x 
#' @rdname Post_x
#' @title Posterior predictive distribution of the count in the treated group
#' @description Density, distribution function, quantile function and random 
#' generation for the posterior predictive distribution of the count in the treated group.
#' @details The posterior predictive distribution of the count in the treated group is a  
#' \code{\link[=PGIBDist]{Poisson-Gamma-Inverse Beta distribution}}.
#' 
#' @param xnew,q vector of non-negative \strong{integer} quantiles 
#' @param x,y counts (integer) in the treated group and control group of the observed experiment
#' @param p vector of probabilities
#' @param a non-negative shape parameter of the Gamma prior distribution on the rate \eqn{\mu}
#' @param c,d non-negative shape parameters of the prior distribution on \eqn{\phi} 
#' @param S,Snew sample sizes of the treated group in the observed experiment and the
#' predicted experiment
#' @param n number of observations to be simulated
#' @param ... arguments passed to \code{\link{summary_PGIB}}
#' 
#' @return \code{dpost_x} gives the density, \code{ppost_x} the distribution function, 
#' \code{qpost_x} the quantile function, \code{rpost_x} samples from the distribution, 
#' and \code{spost_x} gives a summary of the distribution.
#' 
#' @note \code{Post_x} is a generic name for the functions documented. 
#' 
#' @examples 
#' barplot(dpost_x(0:10, 10, 2, 3, 4, 5, 3, 10))
#' qpost_x(0.5, 10, 2, 3, 4, 5, 3, 10)
#' ppost_x(4, 10, 2, 3, 4, 5, 3, 10)
#' @importFrom stats setNames
NULL
#'
#' @rdname Post_x
#' @export 
dpost_x <- function(xnew, Snew, a=0.5, c=0.5, d=0, x, y, S){
  a.post <- a+x+y
  c.post <- c+x
  d.post <- d+a+y
  return( setNames(dPGIB(xnew,a.post,d.post,c.post,S/Snew), paste("xnew=",xnew,sep="")) )
}
#'
#' @rdname Post_x
#' @export 
ppost_x <- function(q, Snew, a=0.5, c=0.5, d=0, x, y, S){
  a.post <- a+x+y
  c.post <- c+x
  d.post <- d+a+y
  return( setNames(pPGIB(q,a.post,d.post,c.post,S/Snew), paste("xnew\u2264",q,sep="")) )
}
#'
#' @rdname Post_x
#' @export 
qpost_x <- function(p, Snew, a=0.5, c=0.5, d=0, x, y, S){
  a.post <- a+x+y
  c.post <- c+x
  d.post <- d+a+y
  return( qPGIB(p,a.post,d.post,c.post,S/Snew) ) 
}
#'
#' @rdname Post_x
#' @export 
rpost_x <- function(n, Snew, a=0.5, c=0.5, d=0, x, y, S){
  a.post <- a+x+y
  c.post <- c+x
  d.post <- d+a+y
  return( rPGIB(n,a.post,d.post,c.post,S/Snew) )
}
#'
#' @rdname Post_x
#' @export 
spost_x <- function(Snew, a=0.5, c=0.5, d=0, x, y, S, ...){
  a.post <- a+x+y
  c.post <- c+x
  d.post <- d+a+y
  return( summary_PGIB(a.post, d.post, c.post, S/Snew, ...) )
}


#' @name Post_y 
#' @rdname Post_y
#' @title Posterior predictive distribution of the count in the control group
#' @description Density, distribution function, quantile function and random 
#' generation for the posterior predictive distribution of the count in the control group.
#' @details The posterior predictive distribution of the count in the treated group is a  
#' \code{\link[=PGIBDist]{Poisson-Gamma-Inverse Beta distribution}}.
#' 
#' @param ynew,q vector of non-negative \strong{integer} quantiles 
#' @param x,y counts (integer) in the treated group and control group of the observed experiment
#' @param p vector of probabilities
#' @param a,b non-negative shape parameter and rate parameter of the Gamma prior distribution on the rate \eqn{\mu}
#' @param c,d non-negative shape parameters of the prior distribution on \eqn{\phi} 
#' @param T,Tnew sample sizes of the control group in the observed experiment and the
#' predicted experiment
#' @param n number of observations to be simulated
#' @param ... arguments passed to \code{\link{summary_PGIB}}
#' 
#' @return \code{dpost_y} gives the density, \code{ppost_y} the distribution function, 
#' \code{qpost_y} the quantile function, \code{rpost_y} samples from the distribution, 
#' and \code{spost_y} gives a summary of the distribution.
#' 
#' @note \code{Post_y} is a generic name for the functions documented. 
#' 
#' @examples 
#' barplot(dpost_y(0:10, 10, 2, 7, 3, 4, 5, 3, 10))
#' spost_y(10, 2, 7, 3, 4, 5, 3, 10, output="pandoc")
#' @importFrom stats setNames
NULL
#'
#' @rdname Post_y
#' @export 
dpost_y <- function(ynew, Tnew, a=0.5, b=0, c=0.5, d=0, x, y, T){
  a.post <- a+x+y
  c.post <- c+x
  d.post <- d+a+y
  return( setNames(dPGIB(ynew,a.post,c.post,d.post,(b+T)/Tnew), paste("ynew=",ynew,sep="")) )
}
#'
#' @rdname Post_y
#' @export 
ppost_y <- function(q, Tnew, a=0.5, b=0, c=0.5, d=0, x, y, T){
  a.post <- a+x+y
  c.post <- c+x
  d.post <- d+a+y
  return( setNames(pPGIB(q,a.post,c.post,d.post,(b+T)/Tnew), paste("ynew\u2264",q,sep="")) )
}
#'
#' @rdname Post_y
#' @export 
qpost_y <- function(p, Tnew, a=0.5, b=0, c=0.5, d=0, x, y, T){
  a.post <- a+x+y
  c.post <- c+x
  d.post <- d+a+y
  return( qPGIB(p,a.post,c.post,d.post,(b+T)/Tnew) ) 
}
#'
#' @rdname Post_y
#' @export 
rpost_y <- function(n, Tnew, a=0.5, b=0, c=0.5, d=0, x, y, T){
  a.post <- a+x+y
  c.post <- c+x
  d.post <- d+a+y
  return( rPGIB(n,a.post,c.post,d.post,(b+T)/Tnew) )
}
#'
#' @rdname Post_y
#' @export 
spost_y <- function(Tnew, a=0.5, b=0, c=0.5, d=0, x, y, T, ...){
  a.post <- a+x+y
  c.post <- c+x
  d.post <- d+a+y
  return( summary_PGIB(a.post,c.post,d.post,(b+T)/Tnew, ...) )
}
