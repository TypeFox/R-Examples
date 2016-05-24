#' Summary of a Gamma distribution
#' 
#' Mode, mean, variance, and quartiles for a Gamma distribution 
#' with shape parameter \code{a} and rate parameter \code{b}.
#'
#' @param a,b Shape and rate parameters.
#' @param output \code{"list"} to return a list, \code{"pandoc"} to print a table
#' @param ... arguments passed to \code{\link[=pander]{pander.data.frame}}
#' @examples
#' summary_gamma(a=2, b=4, output="pandoc", style="rmarkdown")
#' @importFrom pander pander
#' @importFrom stats qgamma
#' @export
summary_gamma <- function(a, b, output="list", ...){
  out <- list(mode=ifelse(a>1, (a-1)/b, 0),  
       mean=a/b,  
       sd=sqrt(a)/b,  
       Q1=qgamma(0.25,a,b),
       Q2=qgamma(0.5,a,b),
       Q3=qgamma(0.75,a,b)
  )
  if(output=="pandoc"){
    pander(data.frame(out), ...)
    return(invisible())
  }else{
    return(out)
  }
}

#' Summary of a Negative Binomial distribution
#' 
#' Mode, mean, variance, and quartiles for a Negative Binomial distribution 
#' with shape parameter \code{size} and probability parameter \code{prob}.
#'
#' @param size,prob parameters of the negative binomial distribution (as for \code{\link[stats]{NegBinomial}})
#' @param output \code{"list"} to return a list, \code{"pandoc"} to print a table
#' @param ... arguments passed to \code{\link[=pander]{pander.data.frame}}
#'
#' @examples
#' summary_nbinom(size=2, prob=0.4, output="pandoc", style="rmarkdown")
#' @importFrom pander pander
#' @importFrom stats qnbinom
#' @export
summary_nbinom <- function(size, prob, output="list", ...){
  out <- list(mode=ifelse(size>1, floor((1-prob)*(size-1)/prob), 0),  
              mean=size*(1-prob)/prob,  
              sd=sqrt(size*(1-prob))/prob,  
              Q1=qnbinom(0.25, size, prob),
              Q2=qnbinom(0.5, size, prob),
              Q3=qnbinom(0.75, size, prob)
  )
  if(output=="pandoc"){
    pander(data.frame(out), ...)
    return(invisible())
  }else{
    return(out)
  }
}


#' @name Beta2Dist 
#' @rdname Beta2Dist
#' @title Beta distribution of the second kind
#' @description Density, distribution function, quantile function and random 
#' generation for the Beta distribution of the second kind with shape parameters 
#' \code{c} and \code{d} and scale parameter \code{scale}. 
#' @details The Beta distribution of the second kind with shape parameters 
#' \eqn{c>0} and \eqn{d>0} and scale parameter \eqn{k>0} is the distribution of 
#' \eqn{k*(U/(1-U))} where \eqn{U} is a random variable following the Beta distribution 
#' with shape parameters  \eqn{c} and \eqn{d}. 
#' 
#' @param x,q vector of quantiles 
#' @param p vector of probabilities
#' @param c,d non-negative shape parameters
#' @param scale non-negative scale parameter
#' @param n number of observations to be simulated
#' @param log logical; if true, returns the logarithm of the result
#' @param output type of the \code{summary_beta2} output: \code{"list"} to return a list, \code{"pandoc"} to print a table
#' @param ... other arguments passed to \code{\link{FDist}}
#' 
#' @return \code{dbeta2} gives the density, \code{pbeta2} the distribution function, \code{qbeta2} the quantile function, and \code{rbeta2} generates random observations, 
#' and \code{summary_beta2} returns a summary of the distribution.
#' 
#' @note \code{Beta2Dist} is a generic name for the functions documented. 
#' 
#' @importFrom pander pander
#' @importFrom stats df pf qf rf
#' @examples 
#' curve(dbeta2(x, 3, 10, scale=2), from=0, to=3)
#' u <- rbeta(1e5, 3, 10)
#' lines(density(2*u/(1-u)), col="blue", lty="dashed")
#' summary_beta2(3,10,2)
#' 
NULL
#'
#' @rdname Beta2Dist
#' @export 
dbeta2 <- function(x, c, d, scale, log=FALSE, ...){
  k <- d/c/scale
  switch(as.character(log),
         "FALSE" = k*df(k*x, df1=2*c, df2=2*d, log=FALSE, ...), 
         "TRUE" = log(k)+df(k*x, df1=2*c, df2=2*d, log=TRUE, ...)
  )
}
#' 
#' @rdname Beta2Dist 
#' @export 
pbeta2 <- function(q, c, d, scale, ...){
  pf(d/c/scale*q, df1=2*c, df2=2*d, ...)
}
#'
#' @rdname Beta2Dist 
#' @export
qbeta2 <- function(p, c, d, scale, ...){
  scale*c/d*qf(p, df1=2*c, df2=2*d, ...)
}
#'
#' @rdname Beta2Dist 
#' @export 
rbeta2 <- function(n, c, d, scale){
  scale*c/d*rf(n, df1=2*c, df2=2*d)
}
#' @rdname Beta2Dist
#' @export
summary_beta2 <- function(c, d, scale, output="list", ...){
  out <- list(mode=ifelse(c>1, scale*(c-1)/(d+1), 0),  
       mean=ifelse(d>1, scale*c/(d-1), Inf),  
       sd=ifelse(d>2, sqrt(scale^2*c*(c+d-1)/((d-1)^2*(d-2))), Inf),  
       Q1=qbeta2(0.25,c,d,scale),
       Q2=qbeta2(0.5,c,d,scale),
       Q3=qbeta2(0.75,c,d,scale)
  )
  if(output=="pandoc"){
    pander(data.frame(out), ...)
    return(invisible())
  }else{
    return(out)
  }
}


#' @name GIBDist 
#' @rdname GIBDist
#' @title Gamma-Inverse Beta distribution
#' @description Density and random  generation for the Gamma-Inverse Beta distribution 
#' with shape parameters \code{a}, \code{alpha}, \code{beta} and rate parameter \code{rho}. 
#' @details This is the mixture distribution obtained by sampling a value \eqn{b} from a Beta distribution 
#' with shape parameters \eqn{\beta}, \eqn{\alpha} 
#' and then sampling a Gamma distribution with shape \eqn{a} and rate \eqn{\rho/b}.
#' 
#' @param x vector of quantiles
#' @param a non-negative shape parameter of the Gamma distribution
#' @param alpha,beta non-negative shape parameters of the mixing Beta distribution
#' @param rho rate parameter 
#' @param n number of observations to be simulated
#' @param output type of the \code{summary_GIB} output: \code{"list"} to return a list, \code{"pandoc"} to print a table
#' @param ... arguments passed to \code{\link[=pander]{pander.data.frame}}
#' 
#' @return \code{dGIB} gives the density, \code{rGIB} samples from the distribution, 
#' and \code{summary_GIB} returns a summary of the distribution.
#' 
#' @note \code{GIBDist} is a generic name for the functions documented. 
#' 
#' @importFrom gsl lnpoch lngamma hyperg_U
#' @importFrom pander pander
#' @importFrom stats rgamma rbeta
#' @examples
#' curve(dGIB(x, 3, 4, 2, 2.5), from=0, to=3)
#' summary_GIB(3, 4, 2, 2.5, output="pandoc", style="grid")
#' 
NULL
#'
#' @rdname GIBDist
#' @export
dGIB <- function(x,a,alpha,beta,rho){
  exp(lnpoch(beta,alpha)-lngamma(a))*rho^a*x^(a-1)*exp(-rho*x)*hyperg_U(alpha,a-beta+1,rho*x)
}
#'
#' @rdname GIBDist
#' @export
rGIB <- function(n,a,alpha,beta,rho){
  rbeta(n,beta,alpha)*rgamma(n, a, rho)
}
#'
#' @rdname GIBDist
#' @export
summary_GIB <- function(a, alpha, beta, rho, output="list", ...){
  out <- list(mean=a*beta/(alpha+beta)/rho,
       sd = sqrt(a*(1+a)*beta*(beta+1)/(alpha+beta)/(alpha+beta+1)/rho^2 - (a*beta/(alpha+beta)/rho)^2)
       )
  if(output=="pandoc"){
    pander(data.frame(out), ...)
    return(invisible())
  }else{
    return(out)
  }
}


#' @name PGIBDist 
#' @rdname PGIBDist
#' @title Poisson-Gamma-Inverse Beta distribution
#' @description Density and random  generation for the Poisson-Gamma-Inverse Beta distribution 
#' with shape parameters \code{a}, \code{c}, \code{d} and scale parameter \code{rho}. 
#' @details This is the mixture distribution obtained by sampling a value from a 
#' \link[=GIBDist]{Gamma-Inverse Beta distribution} and then sampling from 
#' a Poisson distribution having this value as mean.
#' 
#' @param x,q vector of \strong{integer} quantiles
#' @param p vector of probabilities
#' @param a non-negative shape parameter of the Gamma distribution
#' @param alpha,beta non-negative shape parameters of the mixing Beta distribution
#' @param rho hyperrate parameter (rate of the mixing distribution)
#' @param output type of the \code{summary_PGIB} output: \code{"list"} to return a list, \code{"pandoc"} to print a table
#' @param n number of observations to be simulated
#' @param ... arguments passed to \code{\link[=pander]{pander.data.frame}}
#' 
#' @return \code{dPGIB} gives the density, \code{rPGIB} samples from the distribution, 
#' and \code{summary_PGIB} gives a summary of the distribution.
#' 
#' @note \code{PGIBDist} is a generic name for the functions documented. 
#' 
#' @importFrom gsl lnpoch lngamma lnfact
#' @importFrom pander pander
#' @importFrom stats rpois
#' @examples
#' barplot(dPGIB(0:5, a=13, alpha=4, beta=2, rho=2.5), names=0:5)
#' summary_PGIB(13, 4, 2, 2.5)
NULL
#'
#' @rdname PGIBDist
#' @export
dPGIB <- function(x,a,alpha,beta,rho){
    exp(lnpoch(a,x)+lnpoch(beta,x)-lnpoch(alpha+beta,x)-lnfact(x)-x*log(rho))*
      Gauss2F1(beta+x, a+x, alpha+beta+x, -1/rho)
}
#'
#' @rdname PGIBDist
#' @export
pPGIB <- function(q, a, alpha, beta, rho){
  return( vapply(q, FUN=function(x) sum(dPGIB(0:x, a, alpha, beta, rho)), 
                 FUN.VALUE=numeric(1)) )
}
#'
#' @rdname PGIBDist
#' @export
qPGIB <- function(p, a, alpha, beta, rho){
  return( vapply(p, FUN=function(x) icdf(dPGIB, x, a=a, alpha=alpha, beta=beta, rho=rho),
                 FUN.VALUE=numeric(1)) )
}
#'
#' @rdname PGIBDist
#' @export
rPGIB <- function(n, a, alpha, beta, rho){
  return( rpois(n, rGIB(n, a, alpha, beta, rho)) )
}
#'
#' @rdname PGIBDist
#' @export
summary_PGIB <- function(a, alpha, beta, rho, output="list", ...){
  out <- list(mean=a*beta/(alpha+beta)/rho,
              sd = sqrt( a*beta/rho*( 1/(alpha+beta) + (alpha*(a+beta+1)+beta*(beta+1))/(alpha+beta)^2/(alpha+beta+1)/rho )),
              Q1 = qPGIB(0.25, a, alpha, beta, rho),
              Q2 = qPGIB(0.5, a, alpha, beta, rho),
              Q3 = qPGIB(0.75, a, alpha, beta, rho)
  )
  if(output=="pandoc"){
    pander(data.frame(out), ...)
    return(invisible())
  }else{
    return(out)
  }
}

#' @name BNBDist
#' @rdname BNBDist
#' @title Beta-negative binomial distribution
#' @description Density, cumulative function, quantile function and random generation 
#' for the  Beta-negative binomial distribution 
#' with shape parameters \code{a}, \code{c}, \code{d}. 
#' @details This is the mixture distribution obtained by sampling a value \eqn{b} 
#' from a Beta distribution with parameters \eqn{c}, \eqn{d},  
#' then sampling a value \eqn{\lambda} from a Gamma distribution with 
#' shape \eqn{a} and rate \eqn{b/(1-b)}, and
#' then sampling a Poisson distribution with mean \eqn{\lambda}.
#' 
#' @param x,q vector of non-negative integer quantities
#' @param p vector of probabilities
#' @param a,c,d non-negative shape parameters
#' @param n number of observations to be sampled
#' @param output type of the \code{summary_beta_nbinom} output: \code{"list"} to return a list, \code{"pandoc"} to print a table
#' @param ... other arguments passed to \code{\link[SuppDists]{ghyper}}
#' 
#' @return \code{dbeta_nbinom} gives the density, 
#' \code{pbeta_nbinom} the cumulative function, 
#' \code{qbeta_nbinom} the quantile function, 
#' \code{rbeta_nbinom} samples from the distribution, 
#' \code{sbeta_nbinom} and \code{summary_beta_nbinom} give some summaries of the distribution.
#' 
#' @note \code{BNBDist} is a generic name for the functions documented. 
#' 
#' @importFrom SuppDists dghyper pghyper qghyper rghyper sghyper
#' @importFrom pander pander
#' @examples
#' a <- 2 ; c <- 5 ; d <- 30
#' barplot(dbeta_nbinom(0:50, a, c, d), names=0:50)
#' summary_beta_nbinom(a, c, d)
NULL
#'
#' @rdname BNBDist
#' @export
dbeta_nbinom <- function(x, a, c, d, ...){
  dghyper(x, -d, -a, c-1, ...)
}
#
#' @rdname BNBDist
#' @export
pbeta_nbinom <- function(q, a, c, d, ...){
  pghyper(q, -d, -a, c-1, ...)
}
#'
#' @rdname BNBDist
#' @export 
qbeta_nbinom <- function(p, a, c, d, ...){
  qghyper(p, -d, -a, c-1, ...)
}
#'
#' @rdname BNBDist
#' @export 
rbeta_nbinom <- function(n, a, c, d){
  rghyper(n, -d, -a, c-1)
}
#'
#' @rdname BNBDist
#' @export 
sbeta_nbinom <- function(a, c, d){
  sghyper(-d, -a, c-1)
}
#'
#' @rdname BNBDist
#' @export 
summary_beta_nbinom <- function(a, c, d, output="list", ...){
  out <- c(with(sghyper(-d, -a, c-1), 
              list(mean=Mean,
                   mode=Mode,
                   sd=SD)),
              list(
                   Q1=qbeta_nbinom(.25, a, c, d),
                   Q2=qbeta_nbinom(.5, a, c, d),
                   Q3=qbeta_nbinom(.75, a, c, d)
                   )
           )
  if(output=="pandoc"){
    pander(data.frame(out), ...)
    return(invisible())
  }else{
    return(out)
  }
}


#' @name GB2Dist
#' @rdname GB2Dist
#' @title Gamma-Beta2 distribution
#' @description Density and random generation for the  Gamma-Beta2 distribution
#' with shape parameters \code{a}, \code{c}, \code{d} 
#' and rate parameter \code{tau} (scale of the Beta2 distribution). 
#' @details This is the mixture distribution obtained by sampling a value \eqn{y} 
#' from the \link[=Beta2Dist]{Beta2 distribution} with shape parameters \eqn{c}, \eqn{d}, 
#' and scale \eqn{\tau} and 
#' then sampling a value  from the Gamma distribution with 
#' shape \eqn{a} and rate \eqn{y}.
#' The pdf involves 
#' the Kummer confluent hypergeometric function of the second kind. 
#' The cdf involves the generalized hypergeometric function. Its current implementation 
#' does not work when \code{a-d} is an integer, and also fails for many other cases.
#' 
#' @param x,q vector of non-negative quantiles
#' @param p vector of probabilities
#' @param a,c,d non-negative shape parameters
#' @param tau non-negative rate parameter 
#' @param k the order of the moment
#' @param output type of the \code{summary_GB2} output: \code{"list"} to return a list, \code{"pandoc"} to print a table
#' @param ... arguemnts passed to \code{\link[=hypergeo]{genhypergeo}} function
#' @param n number of observations to be sampled
#' 
#' @return \code{dGB2} gives the density, \code{pGB2} the cumulative function, 
#' \code{rGB2} samples from the distribution, and \code{summary_GB2} gives a summary 
#' of the distribution.
#' 
#' @note \code{GB2Dist} is a generic name for the functions documented. 
#' 
#' @importFrom gsl lnpoch lnbeta hyperg_U
#' @importFrom pander pander
#' @importFrom stats rgamma
#' 
#' @examples
#' a <- 2 ; c <- 4 ; d <- 3; tau <- 1.67
#' sims <- rGB2(1e6, a, c, d, tau)
#' mean(sims); moment_GB2(1,a,c,d,tau)
#' mean(sims^2); moment_GB2(2,a,c,d,tau)
#' summary_GB2(a,c,d,tau)
NULL
#'
#' @rdname GB2Dist
#' @export
dGB2 <- function(x,a,c,d,tau){
  return(
    ifelse(d>1 & x<.Machine$double.eps, 0, 
           tau^a*exp(lnpoch(a,c)-lnbeta(d,c)+(a-1)*log(x)+log(hyperg_U(a+c,a-d+1,tau*x)))
    )
  )
}
#'
#' @rdname GB2Dist
#' @export
pGB2 <- function(q, a, c, d, tau, ...){
  return(  exp(a*log(tau) + lnpoch(a,c) - lnbeta(d,c))*
             (sign(d-a)*exp(a*log(q)+lngamma(d-a)-lngamma(c+d)-log(a))*
                genhypergeo(U=c(a,a+c), L=c(1+a,1+a-d), tau*q, ...) +
                sign(a-d)*exp((d-a)*log(tau)+d*log(q)+lngamma(a-d)-lngamma(a+c)-log(d))*
                genhypergeo(U=c(d,c+d), L=c(1+d,1+d-a), tau*q, ...)
             )
  )
}
#'
#' @rdname GB2Dist
#' @export
qGB2 <- function(p, a, c, d, tau){
  return( icdf(dGB2, p, a=a, c=c, d=d, tau=tau) )
}
#
#' @rdname GB2Dist
#' @export
rGB2 <- function(n, a, c, d, tau){
  return( rgamma(n, a, rbeta2(n, c, d, tau)) )
}
#'
#' @rdname GB2Dist
#' @export
moment_GB2 <- function(k,a,c,d,tau){
  return(
    exp(-k*log(tau)+lnpoch(a,k)+lnpoch(d,k)-lnpoch(c-k,k))
  )
}
#'
#' @rdname GB2Dist
#' @export
summary_GB2 <- function(a, c, d, tau, output="list", ...){
  out <- list(mean=ifelse(c>1, moment_GB2(1, a, c, d, tau), Inf),
              sd=ifelse(c>2, sqrt(moment_GB2(2, a, c, d, tau) - (moment_GB2(1, a, c, d, tau))^2), Inf)
                )
  if(output=="pandoc"){
    pander(data.frame(out), ...)
    return(invisible())
  }else{
    return(out)
  }
}


#' @name PGB2Dist
#' @rdname PGB2Dist
#' @title Poisson-Gamma-Beta2 distribution
#' @description Density and random generation 
#' for the  Poisson-Gamma-Beta2 distribution
#' with shape parameters \code{a}, \code{c}, \code{d} 
#' and hyperrate parameter \code{tau} (scale of the Beta2 distribution). 
#' For \code{tau=1} this is the same as the \link[=BNBDist]{Beta-negative binomial distribution}.
#' @details This is the mixture distribution obtained by sampling a value \eqn{y} 
#' from the \link[=Beta2Dist]{Beta2 distribution} with shape parameters \eqn{c}, \eqn{d}, 
#' and scale \eqn{\tau},   
#' then sampling a value \eqn{\lambda} from the Gamma distribution with 
#' shape \eqn{a} and rate \eqn{y}, and
#' then sampling the Poisson distribution with mean \eqn{\lambda}.
#' 
#' @param x,q vector of non-negative \strong{integer} quantiles
#' @param p vector of probabilities
#' @param a,c,d non-negative shape parameters
#' @param tau non-negative hyperrate parameter 
#' @param output type of the \code{summary_PGB2} output: \code{"list"} to return a list, \code{"pandoc"} to print a table
#' @param n number of observations to be sampled
#' 
#' @return \code{dPGB2} gives the density, \code{pPGB2} the cumulative function, 
#' \code{rPGB2} samples from the distribution, and \code{summary_PGB2} gives 
#' a summary of the distribution. 
#' 
#' @importFrom stats rnbinom
#' 
#' @note \code{PGB2Dist} is a generic name for the functions documented. 
#' 
#' @examples
#' a <- 2 ; c <- 5 ; d <- 30; tau <- 2
#' barplot(dPGB2(0:40, a, c, d, tau), names=0:40)
#' summary_PGB2(a,c,d,tau, output="pandoc")
NULL
#'
#' @rdname PGB2Dist
#' @export
dPGB2 <- function(x,a,c,d,tau){
  return( exp(a*log(tau) + dbeta_nbinom(x,a,c,d,log=TRUE) + log(Gauss2F1(a+c,a+x,a+x+c+d,1-tau))) )
}
#'
#' @rdname PGB2Dist
#' @export
pPGB2 <- function(q, a, c, d, tau){
  return( vapply(q, FUN=function(x)  sum(dPGB2(0:x, a, c, d, tau)), 
                 FUN.VALUE=numeric(1)) )
}
#'
#' @rdname PGB2Dist
#' @export
qPGB2 <- function(p, a, c, d, tau){
  return( vapply(p, FUN=function(x) icdf(dPGB2, x, a=a, c=c, d=d, tau=tau),
                 FUN.VALUE=numeric(1)) )
}
#
#' @rdname PGB2Dist
#' @export
rPGB2 <- function(n, a, c, d, tau){
  psi <- rbeta2(n, c, d, tau)
  return( rnbinom(n, a, psi/(1+psi)) )
}
#'
#' @rdname PGB2Dist
#' @export
summary_PGB2 <- function(a, c, d, tau, output="list"){
  m1 <- moment_GB2(1, a, c, d, tau)
  out <- list(mean = ifelse(c>1, m1, Inf),
              sd = ifelse(c>2, sqrt(m1 + moment_GB2(2, a, c, d, tau) - m1^2)),
              Q1 = qPGB2(0.25, a, c, d, tau),
              Q2 = qPGB2(0.5, a, c, d, tau),
              Q3 = qPGB2(0.75, a, c, d, tau)
  )
  if(output=="pandoc"){
    pander(data.frame(out))
    return(invisible())
  }else{
    return(out)
  }
}

