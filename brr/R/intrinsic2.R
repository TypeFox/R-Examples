#' Second intrinsic discrepancy
#' 
#' Intrinsic discrepancy from \code{phi0} to \code{phi} in the marginal model.
#' 
#' @param phi0 the proxy value of \code{phi}
#' @param phi the true value of the parameter
#' @param a,b, the parameters of the prior Gamma distribution on \eqn{\mu}
#' @param S,T sample sizes
#' @return A number, the intrinsic discrepancy from \code{phi0} to \code{phi}.
#' @export
intrinsic2_discrepancy <- function(phi0, phi, a, b, S, T){
  tmp <- phi
  phi <- pmin(phi,phi0)
  phi0 <- pmax(tmp,phi0)
  bphi <- phi*S/(T+b)
  bphi0 <- phi0*S/(T+b)
  N <- log1p(bphi0) - log1p(bphi) # log((bphi0+1)/(bphi+1))
return( a/b*(T+b)*(N+ifelse(bphi<.Machine$double.eps^5, 0, bphi*(N+log(bphi/bphi0)))) ) 
}


#' @name Intrinsic2Inference
#' @rdname Intrinsic2Inference
#' @title Intrinsic inference on the rates ratio based on the second intrinsic discrepancy.
#' 
#' @param a,b,c,d Prior parameters
#' @param S,T sample sizes
#' @param x,y Observed counts
#' @param phi0 the proxy value of \code{phi}
#' @param beta_range logical, if \code{TRUE} (default), an internal method is used to avoid a possible failure in numerical integration; see the main vignette for details
#' @param nsims number of simulations
#' @param conf credibility level
#' @param phi.star the hypothesized value of \code{phi} 
#' @param alternative alternative hypothesis, "less" for H1: \code{phi0 < phi.star}, 
#' "greater" for  H1: \code{phi0 > phi.star} 
#' @param parameter parameter of interest: relative risk \code{"phi"} or vaccine efficacy \code{"VE"}
#' @param tol accuracy requested 
#' @param ... arguments passed to \code{\link{integrate}}
#' @param otol desired accuracy for optimization
#'
#' @return \code{intrinsic2_phi0} returns the posterior expected loss, 
#' \code{intrinsic2_estimate} returns the intrinsic estimate, 
#' \code{intrinsic2_H0} performs intrinsic hypothesis testing, and 
#' \code{intrinsic2_bounds} returns the intrinsic credibility interval. 
#' 
#' @examples
#' a<-2; b<-10; c<-1/2; d<-0; S<-100; T<-S; x<-0; y<-20
#' intrinsic2_phi0(0.5, x, y, S, T, a, b, c, d)
#' intrinsic2_phi0_sims(0.5, x, y, S, T, a, b, c, d)
#' intrinsic2_estimate(x, y, S, T, a, b, c, d)
#' bounds <- intrinsic2_bounds(x, y, S, T, a, b, c, d, conf=0.95); bounds
#' ppost_phi(bounds[2], a, b, c, d, S, T,  x, y)- ppost_phi(bounds[1], a, b, c, d, S, T, x, y)
#' 
#' @importFrom stats dbeta integrate optimize
NULL
#'
#' @rdname Intrinsic2Inference
#' @export
intrinsic2_phi0 <- function(phi0, x, y,  S, T, a, b, c=0.5, d=0, beta_range=TRUE, tol=1e-8, ...){
  post.c <- x+c
  post.d <- y+a+d
  lambda <- (T+b)/S
  return( vapply(phi0, FUN = function(phi0){ 
    f <- function(u) intrinsic2_discrepancy(phi0, lambda * u/(1-u), a=a, b=b, S=S, T=T)
    range <- if(beta_range) beta_integration_range(post.c, post.d, f, accuracy=tol) else c(0,1)
    integrande <- function(u){
      return( f(u)*dbeta(u, post.c, post.d) )
    }
    I <- integrate(integrande, range[1], range[2], ...)
    return(I$value)
  }, FUN.VALUE=numeric(1)) )
}
#'
#' @rdname Intrinsic2Inference
#' @export
intrinsic2_phi0_sims <- function(phi0, x, y, S, T, a, b, c=0.5, d=0, nsims=1e6){
  post.c <- x+c
  post.d <- y+a+d
  lambda <- (T+b)/S
  sims <- rbeta2(nsims, post.c, post.d, lambda)
  return( vapply(phi0, FUN = function(phi0){ 
    return(mean(intrinsic2_discrepancy(phi0, sims, a=a, b=b, S=S, T=T)))
  }, FUN.VALUE=numeric(1)) )
}
#'
#' @rdname Intrinsic2Inference
#' @export
intrinsic2_estimate <- function(x, y, S, T, a, b, c=0.5, d=0, otol = 1e-8, ...){
  post.cost <- function(u0){
    phi0 <- u0/(1-u0)
    intrinsic2_phi0(phi0, x, y, S, T, a, b, c, d, ...)
  }
  optimize <- optimize(post.cost, c(0, 1), tol=otol)
  u0.min <- optimize$minimum
  estimate <- u0.min/(1-u0.min)
  loss <- optimize$objective
  out <- estimate
  attr(out, "loss") <- loss
  return(out)
}
#'
#' @rdname Intrinsic2Inference
#'  @export 
intrinsic2_H0 <- function(phi.star, alternative, x, y, S, T, a, b, c=0.5, d=0, ...){
  post.c <- x+c
  post.d <- y+a+d
  lambda <- (T+b)/S
  integrande <- function(u){
    intrinsic2_discrepancy(phi.star, lambda * u/(1-u), S=S, T=T, a=a, b=b)*dbeta(u, post.c, post.d)
  }
  psi.star <- phi.star/lambda
  u.star <- psi.star/(1+psi.star)
  bounds <- switch(alternative, less=c(0,u.star), greater=c(u.star, 1))
  value <- integrate(integrande, bounds[1], bounds[2], ...)$value
  return(value)
}
#' 
#' @rdname Intrinsic2Inference 
#' @export
intrinsic2_bounds <- function(x, y, S, T, a, b, c=0.5, d=0, conf=.95, parameter="phi", otol = 1e-08, ...){
  post.cost <- function(phi0){
    intrinsic2_phi0(phi0, x, y, S, T, a, b, c, d, ...)
  }
  post.icdf <- function(p){
    qpost_phi(p, a=a, b=b, c=c, d=d, S=S, T=T, x=x, y=y)
  }
  conf <- min(conf, 1 - conf)
  f <- function(p, post.icdf, conf){
    u.phi <- post.icdf(1 - conf + p) 
    l.phi <- post.icdf(p)
    (post.cost(u.phi)-post.cost(l.phi))^2
  }
  minimize <- optimize(f, c(0, conf), post.icdf = post.icdf, 
                       conf = conf, tol=otol)$minimum
  out <- switch(parameter, 
                phi=c(post.icdf(minimize), post.icdf(1 - conf + minimize)),
                VE = sort(1-c(post.icdf(minimize), post.icdf(1 - conf + minimize))))
  out
}



