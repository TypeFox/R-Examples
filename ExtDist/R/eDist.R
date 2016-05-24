#' @title S3 methods for manipulating eDist objects.
#' @description S3 methods for manipulating eDist objects
#' @rdname eDist
#' @name eDist
#' @aliases logLik.eDist
#' @aliases AIC.eDist
#' @aliases AICc.eDist
#' @aliases MDL.eDist 
#' @aliases vcov.eDist
#' @aliases print.eDist
#' @aliases plot.eDist
#'
#' @note The MDL only works for parameter estimation by numerical maximum likelihood.
#'
#' @param object, x An object of class eDist, usually the output of a parameter estimation function.
#' @param x, A list to be returned as class eDist.
#' @param ... Additional parameters
#' @param k numeric, The penalty per parameter to be used; the default k = 2 is the classical AIC.
#' @param corr logical; should vcov() return correlation matrix (instead of variance-covariance matrix).
#' @param plot logical; if TRUE histogram, P-P and Q-Q plot of the distribution returned else only parameter estimation is returned.
#' @method logLik eDist
#' @method AIC eDist
#' @method AICc eDist
#' @method MDL eDist
#' @method print eDist
#' @method vcov eDist
#' @method plot eDist
#'
#' @author A. Jonathan R. Godfrey, Sarah Pirikahu, and Haizhen Wu.
#'
#' @examples
#' X <- rnorm(20)
#' est.par <- eNormal(X, method ="numerical.MLE")
#' logLik(est.par)
#' AIC(est.par)
#' AICc(est.par)
#' BIC(est.par)
#' MDL(est.par)
#' vcov(est.par)
#' vcov(est.par,corr=TRUE)
#' print(est.par)
#' plot(est.par)

#' @rdname eDist
#' @export
logLik.eDist <- function(object,...){
  lFoo <- get(paste0("l", attributes(object)$distname))
  ll <- lFoo(attributes(object)$ob, w=attributes(object)$weights, params=object)
  return(ll)
}

#' @rdname eDist
#' @export
AIC.eDist <- function(object,..., k=2){
  lFoo <- get(paste0("l", attributes(object)$distname))
  ll <- lFoo(attributes(object)$ob, w=attributes(object)$weights, params=object)
  p <- length(object)
  n <- length(attributes(object)$ob)
  AIC <- k*(-ll + p)
}


#' @rdname eDist
#' @export
AICc <- function(object) UseMethod("AICc")

#' @rdname eDist
#' @export
AICc.eDist <- function(object,...){
  lFoo <- get(paste0("l", attributes(object)$distname))
  ll <- lFoo(attributes(object)$ob, w=attributes(object)$weights, params=object)
  p <- length(object) #number of paramters
  n <- length(attributes(object)$ob)
  AICc <- 2*(-ll + p + p*(p+1)/(n-p-1))
  return(AICc)
}

#' @rdname eDist
#' @export
vcov.eDist <- function(object,..., corr=FALSE){
  vcov = solve(attributes(object)$nll.hessian)
  cor = cov2cor(vcov)
  if(corr){return(cor)} else {return(vcov)}
}

#' @rdname eDist
#' @export
BIC <- function(object) UseMethod("BIC")

#' @rdname eDist
#' @export 
BIC.eDist <- function(object,...){
  lFoo <- get(paste0("l", attributes(object)$distname))
  ll <- lFoo(attributes(object)$ob, w=attributes(object)$weights, params=object)
  n <- length(attributes(object)$ob)
  BIC <- 2*(-ll) + length(object) *log(n) 
  return(BIC)
}

#' @references Myung, I. (2000). The Importance of Complexity in Model Selection. 
#' Journal of mathematical psychology, 44(1), 190-204.

#' @rdname eDist
#' @export
MDL <- function(object) UseMethod("MDL")

#' @rdname eDist
#' @export
MDL.eDist <- function(object,...){
  lFoo <- get(paste0("l", attributes(object)$distname))
  ll <- lFoo(attributes(object)$ob, w=attributes(object)$weights, params=object)  
  nll.hessian <- attributes(object)$nll.hessian
  MDL = ifelse(all(is.finite(nll.hessian)),
               -ll + 1/2*log(abs(det(attributes(object)$nll.hessian))),
               NA)
  return(MDL)
}


#' @rdname eDist
#' @export
print.eDist <- function(x,...){
  cat("\nParameters for the", attributes(x)$distname, "distribution. \n(found using the ", attributes(x)$method, "method.)\n\n")
  if(any(is.na(attributes(x)$par.s.e))) {
    print(data.frame(Parameter=attributes(x)$par.name, 
                     Type=attributes(x)$par.type, 
                     Estimate=attributes(x)$par.vals), 
          row.names=FALSE )
  } else {
    print(data.frame(Parameter=attributes(x)$par.name, 
                     Type=attributes(x)$par.type, 
                     Estimate=attributes(x)$par.vals,
                     S.E. = attributes(x)$par.s.e ), 
          row.names=FALSE )
  }
  cat("\n\n")
}


#' @rdname eDist
#' @export
plot.eDist <- function(x,...){
  data <- attr(x, "ob")
  Dist <- attr(x, "dist")
  est.par <- x
  l <- min(data)
  u <- max(data)
  d <- u-l
  n <- length(data)
  dDist <- function(x) get(paste0("d",Dist))(x,param = est.par)
  pDist <- function(data) get(paste0("p",Dist))(data,param = est.par)
  qDist <- function(data) get(paste0("q",Dist))(data,param = est.par)
    
    op <- par(mfrow=c(2,2)) 
    PerformanceAnalytics::textplot(capture.output(print(x)), valign = "top")
    
    hist(data, probability=TRUE)
    curve(dDist, from=l, to=u, add=TRUE, col="blue")
    
    plot(qDist((1:n-0.5)/n), sort(data), main="Q-Q Plot", xlim = c(l,u), ylim = c(l,u), 
         xlab="Theoretical Quantiles", ylab="Sample Quantiles")
    abline(0,1)
    
    plot((1:n-0.5)/n, pDist(sort(data)), main="P-P Plot", xlim = c(0,1), ylim = c(0,1),
         xlab="Theoretical Percentile", ylab="Sample Percentile")
    abline(0,1)
    
    par(op)
}

