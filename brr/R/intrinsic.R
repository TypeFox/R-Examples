rho <- function(phi, phi0, S, T){
  coef <- S/T
  gamma <- phi*coef
  gamma0 <- phi0*coef
  G <- ifelse(gamma0<=1 & gamma<.Machine$double.eps, log1p(gamma0), gamma * log(gamma/gamma0) + log((gamma0+1)/(gamma+1))*(gamma+1))
  H <- if(gamma0<.Machine$double.eps) gamma else gamma+1 - exp(gamma0/(gamma0+1)*log(gamma/gamma0))*(gamma0+1)
return( pmin(G,H) )
}

# rho <- function(phi, phi0, S, T){
#   coef <- S/T
#   gamma <- phi*coef
#   gamma0 <- phi0*coef
#   G <- gamma * log(gamma/gamma0) + log((gamma0+1)/(gamma+1))*(gamma+1)
#   H <- gamma+1 - exp(gamma0/(gamma0+1)*log(gamma/gamma0))*(gamma0+1)
#   result <- pmin(G,H)
#   result
# }

#' Intrinsic discrepancy
#' 
#' Intrinsic discrepancy from \code{phi0} to \code{(mu,phi)}.
#' 
#' @param phi0 the proxy value of \code{phi}
#' @param mu,phi the true values of the parameters
#' @param S,T sample sizes
#' @return A number, the intrinsic discrepancy from \code{phi0} to \code{(mu,phi)}.
#' @export
intrinsic_discrepancy <- function(phi0, mu, phi, S, T){
  return( mu * T * rho(phi, phi0, S, T) )
}

#' @name IntrinsicInference
#' @rdname IntrinsicInference
#' @title Intrinsic inference on the rate ratio. 
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
#' @param ... other arguments passed to \code{\link{integrate}}
#' @param otol desired accuracy for optimization
#'
#' @return \code{intrinsic_phi0} returns the posterior expected loss, 
#' \code{intrinsic_estimate} returns the intrinsic estimate, 
#' \code{intrinsic_H0} performs intrinsic hypothesis testing, and 
#' \code{intrinsic_bounds} returns the intrinsic credibility interval. 
#' 
#' @examples
#' a<-0.5; b<-0; c<-1/2; d<-0; S<-100; T<-S; x<-0; y<-20
#' intrinsic_phi0(0.5, x, y, S, T, a, b, c, d)
#' intrinsic_phi0_sims(0.5, x, y, S, T, a, b, c, d)
#' intrinsic_estimate(x, y, S, T, a, b, c, d)
#' bounds <- intrinsic_bounds(x, y, S, T, a, b, c, d, conf=0.95); bounds
#' ppost_phi(bounds[2], a, b, c, d, S, T,  x, y)- ppost_phi(bounds[1], a, b, c, d, S, T, x, y)
#' @importFrom stats dbeta integrate optimize
NULL
#'
#' @rdname IntrinsicInference
#' @export
intrinsic_phi0 <- function(phi0, x, y,  S, T, a=0.5, b=0, c=0.5, d=0, beta_range=TRUE, tol=1e-8, ...){
  post.c <- x+c
  post.d <- y+a+d
  post.a <- x+y+a
  lambda <- (T+b)/S
  K <- post.a*post.d/(post.c+post.d)*T/(T+b)
  value <- vapply(phi0, 
                  FUN = function(phi0){
                    f <- function(u) rho(lambda * u/(1-u), phi0, S, T)
                    range <- if(beta_range) beta_integration_range(post.c, post.d+1, f, accuracy=tol) else c(0,1)
                    integrande <- function(u){
                      return( f(u)*dbeta(u, post.c, post.d+1) )
                    }
                    I <- integrate(integrande, range[1], range[2], ...)
                    return(I$value)
                  }, FUN.VALUE=numeric(1))
  return( K*value )
}
#'
#'@rdname IntrinsicInference
#'@export
intrinsic_phi0_sims <- function(phi0, x, y,  S, T, a=0.5, b=0, c=0.5, d=0, nsims=1e6){
  post.c <- x+c
  post.d <- y+a+d
  post.a <- x+y+a
  lambda <- (T+b)/S
  K <-  post.a*post.d/(post.c+post.d)*T/(T+b)
  sims <- rbeta2(nsims, post.c, post.d+1, lambda)
  return( vapply(phi0, FUN = function(phi0){ 
    return(K*mean(rho(sims, phi0, S=S, T=T)))
  }, FUN.VALUE=numeric(1)) )
}
#'
#' @rdname IntrinsicInference
#' @export
intrinsic_estimate <- function(x, y, S, T, a=0.5, b=0, c=0.5, d=0, otol = 1e-08, ...){
  post.cost <- function(u0){
    phi0 <- u0/(1-u0)
    intrinsic_phi0(phi0, x, y, S, T, a, b, c, d, ...)
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
#' @rdname IntrinsicInference
#'  @export 
intrinsic_H0 <- function(phi.star, alternative, x, y, S, T, a=0.5, b=0, c=0.5, d=0, ...){
  post.c <- x+c
  post.d <- y+a+d
  post.a <- x+y+a
  lambda <- (T+b)/S
  K <- post.a*post.d/(post.c+post.d)*T/(T+b)
  integrande <- function(u){
    phi <- lambda * u/(1-u)
    rho(phi, phi.star, S, T)*dbeta(u, post.c, post.d+1)
  }
  psi.star <- phi.star/lambda
  u.star <- psi.star/(1+psi.star)
  bounds <- switch(alternative, less=c(0,u.star), greater=c(u.star, 1))
  value <- integrate(integrande, bounds[1], bounds[2], ...)$value
  return(K*value)
}
#' 
#' @rdname IntrinsicInference 
#' @export
intrinsic_bounds <- function(x, y, S, T, a=0.5, b=0, c=0.5, d=0, conf=.95, parameter="phi", otol = 1e-08, ...){
  post.cost <- function(phi0){
    intrinsic_phi0(phi0, x, y, S, T, a, b, c, d, ...)
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



