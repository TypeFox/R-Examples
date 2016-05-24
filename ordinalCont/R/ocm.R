#' Ordinal regression for continuous scales
#'
#' Continuous ordinal regression with logit link using the 
#' generalized logistic function as g function. 
#' @param formula a formula expression as for regression models, of the form 
#' response ~ predictors. Only fixed effects are supported. 
#' The model must have an intercept: attempts to remove one will lead to a warning and will be 
#' ignored.
#' @param data  an optional data frame in which to interpret the variables occurring in the 
#' formulas
#' @param weights optional case weights in fitting. Defaults to 1.
#' @param start a vector of initial values for the regression coefficients
#' and \code{M},  \code{B}, \code{T}, (offset, slope and symmetry of the g function)
#' @param link link function, i.e. the type of location-scale distribution assumed for the latent 
#' distribution. The default ``logit'' link gives the proportional odds model and is the only link function currently supported.
#' @param gfun A smooth monotonic function capable of capturing the non-linear nature of the 
#' ordinal measure. It defaults to the generalized logistic function, which is currently the only 
#' possibility.
#' @param method The optimizer used to maximize the likelihood function.
#' @keywords likelihood, log-likelihood, ordinal regression.
#' @details Fits a continuous ordinal regression model, with fixed effects. The g function is the generalized logistic function (see \code{\link{g_glf}}), and the link function is the logit, 
#' implying the standard logistic distribution for the latent variable. Maximum likelihood estimation is performed, using \code{optim {stats}} with a quasi-Newton method (\code{"BFGS"}). 
#' For continuous ordinal mixed modelling, see \code{\link{ocmm}}.
#' 
#' @seealso For continuous ordinal mixed models, see \code{\link{ocmm}}
#' @return an object of type \code{ocm} with the components listed below. Parameter estimates are in \code{coefficients}. 
#' The last 3 elements of \code{coefficients} are the parameters of the g function: 
#' \code{M},  \code{B},  and \code{T}.
#' \item{coefficients}{parameter estimates}
#' \item{vcov}{variance-covariance matrix}
#' \item{df}{estimated degrees of freedom}
#' \item{logLik}{value of the log-likelihood at the estimated optimum}
#' \item{len_beta}{number of fixed-effects parameters of the model}
#' \item{len_gfun}{number of parameters in the g function used in the model}
#' \item{fitted.values}{fitted probabilities}
#' \item{residuals}{residuals on the latent scale}
#' \item{v}{vector of continuous scores}
#' \item{x}{model matrix}
#' \item{sample.size}{sample size (can differ from the number of observations if the weights are different from 1)}
#' \item{nobs}{number of observations}
#' \item{call}{call to fit the model}
#' \item{no.pars}{total number of parameters estimated}
#' \item{data}{data frame used}
#' \item{link}{link function used}
#' \item{gfun}{g function used}
#' \item{formula}{formula used}
#'  @references Manuguerra M, Heller GZ (2010). Ordinal Regression Models for Continuous 
#'  Scales, \emph{The International Journal of Biostatistics}: 6(1), Article 14.
#' @author Maurizio Manuguerra, Gillian Heller
#' @export
#' @examples
#' ANZ0001.ocm <- ANZ0001[ANZ0001$cycleno==0 | ANZ0001$cycleno==5,]
#' ANZ0001.ocm$cycleno[ANZ0001.ocm$cycleno==5] <- 1
#' fit.overall  <- ocm(overall  ~ cycleno + age + bsa + treatment, data=ANZ0001.ocm)
#' fit.phys 	  <- ocm(phys 	  ~ cycleno + age + bsa + treatment, data=ANZ0001.ocm)
#' fit.pain 	  <- ocm(pain 	  ~ cycleno + age + bsa + treatment, data=ANZ0001.ocm)
#' fit.mood 	  <- ocm(mood 	  ~ cycleno + age + bsa + treatment, data=ANZ0001.ocm)
#' fit.nausvom  <- ocm(nausvom  ~ cycleno + age + bsa + treatment, data=ANZ0001.ocm)
#' fit.appetite <- ocm(appetite ~ cycleno + age + bsa + treatment, data=ANZ0001.ocm)
#' summary(fit.overall)
#' summary(fit.phys)
#' summary(fit.pain)
#' summary(fit.mood)
#' summary(fit.nausvom)
#' summary(fit.appetite)
#' par(mfrow=c(2,3))
#' plot(fit.overall, CIs='vcov', R=100)
#' plot(fit.phys, CIs='vcov', R=100)
#' plot(fit.pain, CIs='vcov', R=100)
#' plot(fit.mood, CIs='vcov', R=100)
#' plot(fit.nausvom, CIs='vcov', R=100)
#' plot(fit.appetite, CIs='vcov', R=100)
#' par(mfrow=c(1,1))


ocm <- function(formula, data=NULL, weights, start=NULL, link = c("logit"), 
                gfun = c("glf"), method = c("optim", "ucminf"))
{
  if (any(sapply(attributes(terms(formula))$term.labels,function(x)grepl("|", x, fixed=T)))) 
    stop("Random effects specified. Please call ocmm.")
  if (missing(formula)) 
    stop("Model needs a formula")
  if (attributes(terms(formula))$intercept == 0){
    formula <- update(formula, .~.+1)
    warning("The model must have an intercept and it has been added to the formula.")
  }
  link <- match.arg(link)
  gfun <- match.arg(gfun) 
  method <- match.arg(method)
  if(is.null(data)) data <- model.frame(formula=formula, data=parent.frame(n=1))
  if(missing(weights)) weights <- rep(1, nrow(data))
  keep <- weights > 0
  data <- data[keep,]
  weights <- weights[keep]
    
  mf <- model.frame(formula=formula, data=data)
  x <- model.matrix(attr(mf, "terms"), data=mf)
  lenx <- ncol(x)
  v <- model.response(mf)
  xnames <- dimnames(x)[[2]][2:lenx]
  x <- as.matrix(x[,2:lenx])
  colnames(x) <- xnames
  v <- as.numeric(v)
  if (is.null(start)) {
    beta_start <- set.beta_start(x,v)
    len_beta = length(beta_start)
    names(beta_start) <- xnames[1:len_beta]
    if (gfun == 'glf') {
      gfun_start <- set.glf_start(x,v)
      names(gfun_start) <- c("M", "B", "T")
    }
    start <- c(beta_start, gfun_start)
    len_gfun <- length(gfun_start)
  }
  est <- ocmEst(start, v, x, weights, link, gfun, method)
  coef <- est$coefficients
  beta <- coef[1:len_beta]
  par_g <- coef[(len_beta+1):(len_beta+len_gfun)]
  est$len_beta <- len_beta
  est$len_gfun <- len_gfun
  est$fitted.values <- as.vector(-x%*%beta)
  est$residuals <- g_glf(v, par_g) - est$fitted.values
  est$v <- v
  est$x <- x
  est$sample.size <- nrow(x)
  est$nobs <- sum(weights)
  est$call <- match.call()
  est$no.pars <- length(coef)
  est$data <- data
  est$link <- link
  est$gfun <- gfun
  est$method <- method
  est$formula <- formula
  class(est) <- "ocm"
  est
}  

#' @title Log-likelihood function for the fixed-effects model
#'
#' @details  This function computes minus the log-likelihood function for a fixed-effects model using 
#' the generalized logistic function as g function and the logit link function. It is used internally 
#' to fit the model and should not be of interest of the user.
#' @param par vector of regression coefficients, 
#' and \code{M},  \code{B}, \code{T}, (offset, slope and symmetry of the g function)
#' @param v vector of standardized scores from the continuous ordinal scale
#' @param d.matrix design matrix (fixed effects)
#' @param wts optional case weights
#' @param len_beta length of the regression coefficients vector
#' @keywords likelihood, log-likelihood.
#' @return Minus the log-likelihood at parameter values \code{par} 
#' @author Maurizio Manuguerra, Gillian Heller

negloglik_glf <- function(par, v, d.matrix, wts, len_beta){
  return(-sum(wts * logdensity_glf(par, v, d.matrix, len_beta)))
}

logdensity_glf <- function(par, v, d.matrix, len_beta){
  x <- d.matrix
  beta <- par[1:len_beta]
  par_g <- par[(len_beta+1):(len_beta+3)]
  par_dg <- par[(len_beta+2):(len_beta+3)]
  g <- g_glf(v, par_g)
  dg <- dg_glf(v, par_dg)
  if (any(dg<=0)) return(Inf)
  xb <- x %*% beta
  return(log(dg) + g + xb -2*log(1+exp(g+xb)))
}

#' @import ucminf
ocmEst <- function(start, v, x, weights, link, gfun, method){
  len_beta <- ncol(x)
  if (gfun == "glf") {
    if (link == "logit"){
      if (method == "optim"){
        fit <- optim(par=start,negloglik_glf, v=v, d.matrix=x, wts=weights, len_beta=len_beta, method="BFGS", hessian = T)
      } else if (method == "ucminf") {
        fit <- ucminf(par=start,negloglik_glf, v=v, d.matrix=x, wts=weights, len_beta=len_beta, hessian = 3)
      } else {
        stop("Optimization method not implemented.")
      }
    } else {
      stop("link function not implemented.")
    }
  } else {
    stop("g function not implemented.")
  }
  ## compute QR-decomposition of x
  #qx <- qr(x)
  #Hessian
  H=fit$hessian
  #require(numDeriv)
  #H=hessian(negloglik_glf,fit$par,v=v, d.matrix=x,len_beta=len_beta)
  qrH <- qr(H)
  if(qrH$rank < nrow(H))
    stop("Cannot compute vcov: \nHessian is numerically singular")
  vcov <- solve.qr(qrH)
  
  ## compute (x'x)^(-1) x'y
  coef <- fit$par
  names(coef) <- names(start)
  len_beta = ncol(x)
  beta <- coef[1:len_beta]
  par_g <- coef[(len_beta+1):(len_beta+3)]
  
  ## degrees of freedom and standard deviation of residuals
  df <- nrow(x)-ncol(x)-length(par_g)
  fitted.values <- inv.logit(g_glf(v, par_g) + x%*%beta)
  sigma2 <- sum((v - fitted.values)^2)/df
  
  ## compute sigma^2 * (x'x)^-1
  #vcov <- sigma2 * chol2inv(qx$qr)
  colnames(vcov) <- rownames(vcov) <- c(colnames(x),"M", "B", "T")
  list(coefficients = coef,
       vcov = vcov,
##       sigma = sqrt(sigma2),
       df = df,
       logLik = -fit$value)
}
