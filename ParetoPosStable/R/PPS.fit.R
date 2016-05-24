#' @title Fitting the Pareto Positive Stable (PPS) distribution
#' @description \code{PPS.fit()} returns the fit of a PPS distribution to real data, allowing the scale parameter to be held fixed if desired.
#'   @param x a vector of observations
#'   @param estim.method the method of parameter estimation. It may be "MLE", "iMLE", "OLS", or "LMOM".
#'   @param sigma the value of the scale parameter, if it is known; if the value is \code{NULL}, the parameter is estimated.
#'   @param start a list with the initial values of the parameters for some of the estimation methods.
#'   @param Pareto a logical argument to constrain the PPS fit to a Pareto fit when the value is \code{TRUE}.
#'   @param \ldots other arguments.
#'
#' @details 
#' The maximum likelihood method implemented by the direct optimization of the log-likelihood is given by \code{estim.method = "MLE"}. The numerical algorithm to search the optimum is the ``Nelder-Mead'' method implemented in the \code{optim} function, considering as initial values those given in the \code{start} argument or, if it is missing, those provided by the OLS method.
#'   
#' A different approximation of the maximum likelihood estimates is given by \code{estim method = "iMLE"}; it is an iterative methodology where \code{optimize()} function provides the optimum scale parameter value, while the \code{uniroot()} function solve normal equations for that given scale parameter. 
#'   
#' The regression estimates (\code{"OLS"}) searchs an optimum scale value (in a OLS criterion) by the \code{optimize()} function. Then the rest of the parameters are estimated also by OLS, as appears in Sarabia and Prieto (2009).
#'   
#' In the L-moments method (\code{"LMOM"}) estimates are obtained searching parameters that equal the first three sample and theoretical L-moments by means of the ``Nelder-Mead'' algorithm implemented in \code{optim()}; the initial values are given in the \code{start} argument or, if it is missing, provided by the \code{"iMLE"}. 
#' 
#' @return A \code{PPSfit} Object, a list with
#' \itemize{
#'   \item{estimate}{parameter estimates.}
#'   \item{loglik}{the log-likelihood value.}
#'   \item{n}{the number of observations.}
#'   \item{obs}{the observations.}
#'   \item{obsName}{the name of the variable with the observations.}
#'   \item{estim.method}{the method of parameter estimation.}
#'   }
#' When this last value is \code{"LMOM"} the function also returns details about the convergence of the numerical method involved (\code{convergence} value).
#' 
#' @references 
#' Sarabia, J.M and Prieto, F. (2009). The Pareto-positive stable distribution: A new descriptive model for city size data, \emph{Physica A: Statistical Mechanics and its Applications}, \bold{388}(19), 4179-4191.  
#' Hosking, J. R. M. (1990). L-moments: analysis and estimation of distributions using linear combinations of order statistics. \emph{Journal of the Royal Statistical Society, Series B}, \bold{52}, 105-124.
#'   
#' @seealso \code{\link{coef.PPSfit}}, \code{\link{print.PPSfit}}, \code{\link{plot.PPSfit}}, \code{\link{GoF}}
#' 
#' @examples
#' data(turkey)
#' fit <- PPS.fit(turkey$Pop2000)
#' print(fit)
#' coef(fit)
#' se(fit, k = 100, parallel = FALSE)
#' logLik(fit)
#' par(mfrow=c(2,2))
#' plot(fit)
#' GoF(fit, k = 100, parallel = FALSE)

#' @importFrom stats integrate lm nlm optim optimize uniroot

#' @export
PPS.fit <-
  function(x, estim.method = "MLE", sigma = NULL, start, Pareto = FALSE, ...){
    
    if (Pareto==TRUE) pareto.fit(x, estim.method, sigma, start, ...)
    else{
      
      if (!is.element(estim.method,c("MLE", "OLS", "iMLE", "LMOM"))) stop("The estimation method in method.estim is not OLS, MLE, iMLE nor LMOM")
      
      if (missing(x) || length(x) == 0L || mode(x) != "numeric" || any(x<=0)) stop("'x' must be a non-empty numeric vector of positive values")
      if (any(!is.finite(x))) stop("'x' contains missing or infinite values")
      
      if (missing(start)) start <- NULL
      xName <- paste(deparse(substitute(x), 500), collapse = "\n")
      dots <- names(list(...))
      n <- length(x)
      Fn <- ecdf(x)
      y <- log(-log(1-n*Fn(x)/(n+1)))
      
      
      # OLS ---------------------------------------------------------
      OLS.estimate<-function(x, sigma){
        
        if (is.null(sigma)){
          sce <- function(sigma){
            z <- log(x/sigma) 
            v.ols <- sum((log(z)-mean(log(z)))*(y-mean(y)))/sum((log(z)-mean(log(z)))^2)
            lam.ols <- exp(mean(y)-v.ols*mean(log(z)))
            sce <- sum(lm(y~log(z))$residuals^2)
          }
          sigma <- optimize(f=sce,interval=c(0,min(x)))$minimum
          z <- log(x/sigma) 
          v <- sum((log(z)-mean(log(z)))*(y-mean(y)))/sum((log(z)-mean(log(z)))^2)
          lam <- exp(mean(y) - v*mean(log(z)))
          loglik <- n*(log(lam)+log(v))+(v-1)*sum(log(z))-lam*sum(z^v)-sum(log(x))
          return(structure(list(estimate = list(lambda = lam, sigma = sigma, nu = v), loglik = loglik, n = n, obs = x, obsName = xName, estim.method = "OLS"), class = "PPSfit"))
        }
        
        if (!is.null(sigma)){
          sc <- sigma
          z <- log(x/sc)
          v <- sum((log(z)-mean(log(z)))*(y-mean(y)))/sum((log(z)-mean(log(z)))^2)
          lam <- exp(mean(y)-v*mean(log(z)))
          loglik <- n*(log(lam)+log(v))+(v-1)*sum(log(log(x/sc)))-lam*sum(log(x/sc)^v)-sum(log(x))
          return(structure(list(estimate = list(lambda = lam, nu = v), sigma = sigma, loglik = loglik, n = n, obs = x, obsName = xName, estim.method = "OLS"), class = "PPSfit"))
        }
      }
      
      
      # MLE --------------------------------------------------------
      MLE.estimate <- function(x, sigma, start){
        
        if (!is.null(sigma)){
          sc <- sigma
          z <- log(x/sc)
          v.ols <- sum((log(z)-mean(log(z)))*(y-mean(y)))/sum((log(z)-mean(log(z)))^2)
          lam.ols <- exp(mean(y)-v.ols*mean(log(z)))
          obj <- function(v) (1/v+(1/n)*sum(log(z))-sum(z^v*log(z))/sum(z^v))^2
          optimum <- nlm(p = v.ols, f = obj)
          v<-optimum$estimate
          lam <- 1/(sum(z^v)/n)
          loglik <- n*(log(lam)+log(v))+(v-1)*sum(log(log(x/sc)))-lam*sum(log(x/sc)^v)-sum(log(x))
          return(structure(list(estimate = list(lambda = lam, nu = v), sigma = sigma, loglik = loglik, n = n, obs = x, obsName = xName, estim.method = "MLE", convergence = optimum$code), class = "PPSfit"))
        }
        
        if (is.null(sigma)){
          if (is.null(start)){
            sc <- min(x)/2
            Fn <- ecdf(x)
            y <- log(-log(1-n*Fn(x)/(n+1)))
            z <- log(x/sc)
            v <- sum((log(z)-mean(log(z)))*(y-mean(y)))/sum((log(z)-mean(log(z)))^2)
            lam <- exp(mean(y)-v*mean(log(z)))
            start <- list(lambda = lam, sigma = sc, nu = v)
            start <- start[!is.element(names(start), dots)]
          }
          nm <- names(start)
          if (sum(is.element(nm,c("lambda","sigma","nu")))<3) stop("'start' doesn't specify the names of the PPS parameters as 'lambda', 'sigma' and 'nu'")
          
          logL <- function(p){
            lam <- exp(p[1])
            sc <- (1/(1+exp(-p[2])))*min(x)
            v <- exp(p[3])
            -(n*(log(lam)+log(v))+(v-1)*sum(log(log(x/sc)))-lam*sum(log(x/sc)^v)-sum(log(x)))
          }
          p0 <- c(log(start[[1]]), -log((min(x)-start[[2]])/start[[2]]), log(start[[3]]))
          res <- optim(p0, logL, hessian=FALSE)
          l.est <- exp(res$par[1])
          s.est <- (1/(1+exp(-res$par[2])))*min(x)
          v.est <- exp(res$par[3])
          estimate <- list(lambda = l.est, sigma = s.est, nu = v.est)
          return(structure(list(estimate = estimate, loglik = -res$value, n = n, obs = x, obsName = xName, estim.method = "MLE", convergence = res$convergence), class = "PPSfit"))
        }
      }
      
      
      # iMLE --------------------------------------------------------
      iMLE.estimate<-function(x){
        logver <- function(sigma){
          z <- log(x/sigma)
          obj <- function(v) 1/v+(1/n)*sum(log(z))-sum(z^v*log(z))/sum(z^v)
          cota.inf <- 1/(max(log(z))-mean(log(z)))
          cota.sup <- 1/(mean(log(z))-min(log(z)))
          if (cota.inf*cota.sup>=0){
            cota.inf <- 0
            cota.sup <- 100
          }
          v.mle <- uniroot(obj, lower = cota.inf, upper = cota.sup)$root
          lam.mle <- 1/(sum(z^v.mle)/n)
          return(- (n * (log(lam.mle)+log(v.mle)) + (v.mle-1)*sum(log(log(x/sigma))) - lam.mle*sum(log(x/sigma)^v.mle) - sum(log(x))))
        }
        sigma.mle.lv <- optimize(f=logver,interval=c(0,min(x)))$minimum
        z <- log(x/sigma.mle.lv)
        obj <- function(v) 1/v+(1/n)*sum(log(z))-sum(z^v*log(z))/sum(z^v)
        v.mle.lv <- uniroot(obj, lower = 0, upper = 100)$root
        lam.mle.lv <- 1/(sum(z^v.mle.lv)/n)
        loglik.optimum <- n*(log(lam.mle.lv)+log(v.mle.lv))+(v.mle.lv-1)*sum(log(log(x/sigma.mle.lv)))-lam.mle.lv*sum(log(x/sigma.mle.lv)^v.mle.lv)-sum(log(x))
        return(structure(list(estimate = list(lambda = lam.mle.lv, sigma = sigma.mle.lv, nu = v.mle.lv), loglik = loglik.optimum, n = n, obs = x, obsName = xName, estim.method = "iMLE"), class = "PPSfit"))
      }
      
      
      # LMOM -------------------------------------------------------
      LMOM.estimate <- function(x, start){
        code1<-1
        code2<-1
        code3<-1
        l1 <- function(l,s,v) {
          obj1 <- function(p) s*exp((-(1/l)*log(1-p))^(1/v))
          if((integrate(obj1,0,1,stop.on.error = FALSE)$message=="the integral is probably divergent")==TRUE){
            print("Warning: L1 moment evaluated where it is probably divergent")
            return(code1 <- 0)
          }
          else integrate(obj1,0,1)$value
        }
        l2 <- function(l,s,v) {
          obj2 <- function(p) (s*exp((-(1/l)*log(1-p))^(1/v)))*(2*p-1)
          if((integrate(obj2,0,1,stop.on.error = FALSE)$message=="the integral is probably divergent")==TRUE){
            print("Warning: L2 moment evaluated where it is probably divergent")
            return(code2 <- 0)
          }
          else integrate(obj2,0,1)$value
        }
        l3 <- function(l,s,v) {
          obj3 <- function(p) (s*exp((-(1/l)*log(1-p))^(1/v)))*(6*p^2-6*p+1)
          if((integrate(obj3,0,1,stop.on.error = FALSE)$message=="the integral is probably divergent")==TRUE){
            print("Warning: L3 moment evaluated where it is probably divergent")
            return(code3 <- 0)
          }
          else integrate(obj3,0,1)$value
        }
        if ((code1==0)|(code2==0)|(code3==0)) print("The LMOM method has no solution")
        else{
          obj <- function(p) {
            l<-exp(p[1])
            s<-min(x)/(1+exp(-p[2]))
            v<-exp(p[3])
            (lmom::samlmu(x)[[1]]-l1(l,s,v))^2 + (lmom::samlmu(x)[[2]]-l2(l,s,v))^2 + (lmom::samlmu(x)[[3]]-l3(l,s,v))^2
          }
          if (is.null(start)) pars<-iMLE.estimate(x)$estimate
          else pars <- start
          minimo<-optim(par = c(log(pars[[1]]), -log(min(x)/pars[[2]]-1), log(pars[[3]])), fn = obj)
          lam<-exp(minimo$par[1])
          sigma<-min(x)/(1+exp(-minimo$par[2]))
          v<-exp(minimo$par[3])
          minimum <- minimo$value
          loglik <- (n*(log(lam)+log(v))+(v-1)*sum(log(log(x/sigma)))-lam*sum(log(x/sigma)^v)-sum(log(x)))
          return(structure(list(estimate = list(lambda = lam, sigma = sigma, nu = v), minimo = minimum, loglik = loglik, n = n, obs = x, obsName = xName, estim.method = "LMOM",convergence = minimo$convergence), class = "PPSfit"))
        }
      }
      
      
      # Return ------------------------------------------------------------
      if (estim.method == "MLE") return(MLE.estimate(x, sigma, start))
      if (estim.method == "iMLE") return(iMLE.estimate(x))
      if (estim.method == "OLS") return(OLS.estimate(x, sigma))
      if (estim.method == "LMOM") return(LMOM.estimate(x, start))
      #       if (estim.method == "KSMIN") return(KSMIN.estimate(x, start))
    }
  }
