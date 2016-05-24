#' Ordinal regression for continuous scales with mixed effects 
#'
#' Fits an ordinal continuous mixed model with logit link, using the generalized logistic function as g function.
#' @param formula a formula expression as for regression models, of the form 
#' response ~ predictors. Only mixed-effects models with a single random effect on the intercept are supported. 
#' The model must have an intercept: attempts to remove one will lead to a warning and will be 
#' ignored.
#' @param data  an optional data frame in which to interpret the variables occurring in the formulas
#' @param weights optional case weights in fitting. Defaults to 1.
#' @param start a vector of initial values for the regression coefficients, 
#' \code{M},  \code{B}, \code{T}, (offset, slope and symmetry of the g function) and the standard deviation of the random effect
#' @param link link function, i.e., the type of location-scale distribution assumed for the latent distribution. The default logit link gives the proportional odds model and is the only link function currently supported.
#' @param gfun A smooth monotonic function capable of capturing the non-linear nature of the ordinal measure. It defaults to the generalized logistic function (\code{\link{g_glf}}), which is currently the only possibility.
#' @param method The optimizer used to maximize the likelihood function.
#' @param quad A string indicating the type of quadrature used to integrate over the random effects. Can take values \code{"Laplace"} (Adaptive Gauss-Hermite quadrature using Laplace approximation; the default) or \code{"GH"} (Gauss-Hermite quadrature).
#' @param n_nodes order of Gauss-Hermite rule used (number of nodes) 
#'  @details Fits a continuous ordinal regression model, with fixed and random effects. The g function is the generalized logistic function (see \code{\link{g_glf}}), and the link function is the logit, 
#' implying the standard logistic distribution for the latent variable. Maximum likelihood estimation is performed, using \code{optim {stats}} with a quasi-Newton method (\code{"BFGS"}). 
#' Either adaptive Gauss-Hermite quadrature with the Laplace approximation, or Gauss-Hermite quadrature, is used.  
#' For continuous ordinal  modelling with fixed effects only, see \code{\link{ocm}}.
#' @return An object of type \code{ocmm} with the components listed below. 
#' 
#' \item{coefficients}{parameter estimates. The first \code{len_beta} elements
#' are the estimates of the fixed-effects parameters; the last 4  elements  are the estimates of  
#' the parameters of the g function (\code{M},  \code{B},  and \code{T}) and the standard 
#' deviation of the random effect.}
#' \item{vcov}{variance-covariance matrix, of dimension (\code{len_beta} +4)x(\code{len_beta} +4
#' )}
#' \item{sigma_rnd}{standard deviation of the random effect}
#' \item{df}{estimated degrees of freedom}
#' \item{logLik}{value of the log-likelihood at the estimated optimum}
#' \item{len_beta}{number of fixed-effects parameters of the model}
#' \item{len_gfun}{number of parameters in the g function used in the model}
#' \item{len_rnd}{number of random effects (1 in this version of the package)}
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
#' @keywords likelihood, log-likelihood, ordinal regression.
#' @export
#' @author Maurizio Manuguerra, Gillian Heller
#' @examples
#'\dontrun{
#' fit.overall.rnd  <- ocmm(overall  ~ cycleno + age + bsa + treatment + (1|randno), data=ANZ0001)
#' fit.phys.rnd     <- ocmm(phys 	   ~ cycleno + age + bsa + treatment + (1|randno), data=ANZ0001)
#' fit.pain.rnd 	  <- ocmm(pain 	   ~ cycleno + age + bsa + treatment + (1|randno), data=ANZ0001)
#' fit.mood.rnd 	  <- ocmm(mood 	   ~ cycleno + age + bsa + treatment + (1|randno), data=ANZ0001)
#' fit.nausvom.rnd <- ocmm(nausvom ~ cycleno + age + bsa + treatment + (1|randno), data=ANZ0001)
#' fit.appetite.rnd <- ocmm(appetite ~ cycleno + age + bsa + treatment + (1|randno), data=ANZ0001)
#' summary(fit.overall.rnd)
#' summary(fit.phys.rnd)
#' summary(fit.pain.rnd)
#' summary(fit.mood.rnd)
#' summary(fit.nausvom.rnd)
#' summary(fit.appetite.rnd)
#' }



ocmm <- function(formula, data=NULL, weights, start=NULL, link = c("logit"), gfun = c("glf"), method = c("optim", "ucminf"), quad=c("Laplace","GH"), n_nodes=10)
{
  if (missing(formula)) 
    stop("Model needs a formula")
  if (attributes(terms(formula))$intercept == 0){
    formula <- update(formula, .~.+1)
    warning("The model must have an intercept and it has been added to the formula.")
  }
  terms <- attributes(terms(formula))$term.labels
  i_rnd <- which(sapply(attributes(terms(formula))$term.labels,function(x)grepl("|", x, fixed=T)))
  if (length(i_rnd) > 1) stop("Only one random effect supported. Please respecify your model.")
  if (length(i_rnd) == 0) stop("No random effects specified. Please call ocm.")
  #for (rndterms in terms[i_rnd]){ #not needed as only one rnd effect is supported.
  rndterms <- terms[i_rnd]
  fixterms <- terms[-i_rnd]
  ibar <- regexpr("|",rndterms,fixed=T)
  left <- trim(substr(rndterms,1,ibar-1))
  right <- trim(substr(rndterms,ibar+1, nchar(rndterms)))
  #cat("\nLeft:",as.numeric(left)," - Right:",right,"\n")
  if (as.numeric(left)!=1) stop("Only random effects on the intercept are supported in this version of ordinalCont.")
  if (any(sapply(c(":","*","|"), function(x)grepl(x, right,fixed=T)))) stop("Syntax incorrect or feature not implemented.")
  #cat("\nGoing to stratify by",right,"..\n")
  nf.fix <- paste(fixterms, collapse="+")
  nf.rnd <- paste(right, collapse="+")
  form.complete <- update(formula, as.formula(paste("~",nf.fix,"+",nf.rnd)))
  form.fix <- update(formula, as.formula(paste("~",nf.fix)))
  form.rnd <- update(formula, as.formula(paste("~",nf.rnd)))
  
  link <- match.arg(link)
  gfun <- match.arg(gfun)
  method <- match.arg(method)
  if(is.null(data)) data <- model.frame(formula=form.complete, data=parent.frame(n=1))
  quad = match.arg(quad)
  if(missing(weights)) weights <- rep(1, nrow(data))
  keep <- weights > 0
  data <- data[keep,]
  weights <- weights[keep]
  
  to_stratify <- as.factor(data[,right])
  iclusters <- lapply(levels(to_stratify),function(x)which(to_stratify==x))
  
  mf <- model.frame(formula=formula, data=data)
  x.complete <- model.matrix(attr(mf, "terms"), data=mf)
  
  mf.fix <- model.frame(formula=form.fix, data=data)
  mf.rnd <- model.frame(formula=form.rnd, data=data)
  x <- model.matrix(attr(mf.fix, "terms"), data=mf.fix)
  z <- model.matrix(attr(mf.rnd, "terms"), data=mf.rnd)
  xnames <- dimnames(x)[[2]][-1]
  x <- as.matrix(x[,-1]) # 1 for the intercept
  colnames(x) <- xnames
  z <- as.matrix(z)[,-1] # 1 for the intercept
  v <- model.response(mf)
  v <- as.numeric(v)
  if (is.null(start)) {
    beta_start <- set.beta_start(x,v)
    len_beta = length(beta_start)
    names(beta_start) <- xnames[1:len_beta] 
    if (gfun == 'glf') {
      gfun_start <- set.glf_start(x,v)
      names(gfun_start) <- c("M", "B", "T")
    }
    sigma_rnd_eff = 1
    names(sigma_rnd_eff) = paste("Std.Dev.", right, "group")
    start <- c(beta_start, gfun_start, sigma_rnd_eff) #1 is the variance of the single rnd effect
    len_gfun <- length(gfun_start)
  }
  est <- ocmmEst(start, v, x, z, weights, link, gfun, method, rnd=right, n_nodes=n_nodes, quad=quad, iclusters)
  coef <- est$coefficients
  beta <- coef[1:len_beta]
  par_g <- coef[(len_beta+1):(len_beta+len_gfun)]
  est$len_beta <- len_beta
  est$len_gfun <- len_gfun
  est$len_rnd <- 1
  est$rnd <- right
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
  class(est) <- "ocmm"
  est
  
}



#' @title Log-likelihood function for the mixed-effects model
#'@description Log-likelihood function for the mixed-effects continuous ordinal model, using the generalized logistic 
#' function as g function and the logit link function
#' @details This function computes minus the log-likelihood function for a mixed-effects model using the 
#' generalized logistic function as g function and the logit link function. It is used internally to fit the model and should not be of interest of the user.
#' @param par vector of regression coefficients, 
#'\code{M},  \code{B}, \code{T}, (offset, slope and symmetry of the g function)
#' and the standard deviation of the random effect
#' @param v vector of standardized scores from the continuous ordinal scale
#' @param d.matrix design matrix (fixed effects)
#' @param rnd.matrix random term model matrix
#' @param wts optional case weights
#' @param len_beta length of the regression coefficients vector
#' @param rnd character vector listing the random terms
#' @param quad  string indicating the type of quadrature used to integrate over the random 
#' effects. Can take values \code{"Laplace"} (Adaptive Gauss-Hermite quadrature using Laplace 
#' approximation; the default) or \code{"GH"} (Gauss-Hermite quadrature).
#' @param n_nodes order of Gauss-Hermite rule used (number of nodes)
#' @param iclusters  list containing the row numbers of the design matrix relative to each level 
#' of the factor over which random effects are computed
#' @keywords likelihood, log-likelihood.
#' @author Maurizio Manuguerra, Gillian Heller
#' @return Minus the log-likelihood at parameter values \code{par} 

negloglik_glf_rnd <- function(par, v, d.matrix, rnd.matrix, wts, len_beta, rnd, n_nodes, quad, iclusters)
  { 
  #b is the random effect
  densities_by_cluster <- sapply(iclusters, negloglik_glf_rnd2, par=par, v=v, d.matrix=d.matrix, rnd.matrix=rnd.matrix, wts=wts, len_beta=len_beta, n_nodes=n_nodes, quad=quad)
  loglik = -sum(log(densities_by_cluster))
  return(loglik)
}

#' @import fastGHQuad
negloglik_glf_rnd2 <- function(indices, par, v, d.matrix, rnd.matrix, wts, len_beta, n_nodes, quad){
  rule <- gaussHermiteData(n_nodes)
  x <- array(d.matrix[indices,],dim=c(length(indices), len_beta))
  y <- v[indices]
  w <- wts[indices]
  beta <- par[1:len_beta]
  par_g <- par[(len_beta+1):(len_beta+3)]
  par_dg <- par[(len_beta+2):(len_beta+3)]
  sigma_rnd <- par[len_beta+4]
  g <- g_glf(y, par_g)
  dg <- dg_glf(y, par_dg)
  if (any(dg<=0)) return(Inf)
  xb <- as.numeric(x %*% beta)
  if (quad=="Laplace")
    return(aghQuad(g=density_glf_Laplace, muHat=0, sigmaHat=1, rule=rule, all_pars=cbind(g,dg,xb), sigma_rnd=sigma_rnd, w=w))
  else if (quad=="GH")
    return(ghQuad(f=density_glf_GH, rule=rule, all_pars=cbind(g,dg,xb), sigma_rnd=sigma_rnd, w=w))
}

density_glf_Laplace <- function(b, all_pars, sigma_rnd, w){
  lik_cluster <- apply(all_pars, 1, function(x, b, sigma_rnd){x[2]*exp(x[1]+x[3]+sqrt(2)*sigma_rnd*b)/(1+exp(x[1]+x[3]+sqrt(2)*sigma_rnd*b))^2/sqrt(pi)*exp(-b^2)},b=b,sigma_rnd=sigma_rnd)
  lik_cluster <- apply(lik_cluster,1,function(x,w)return(prod(x^w)),w=w)
  return(lik_cluster)
}

density_glf_GH <- function(b, all_pars, sigma_rnd, w){
  lik_cluster <- apply(all_pars, 1, function(x, b, sigma_rnd){x[2]*exp(x[1]+x[3]+sqrt(2)*sigma_rnd*b)/(1+exp(x[1]+x[3]+sqrt(2)*sigma_rnd*b))^2/sqrt(pi)},b=b,sigma_rnd=sigma_rnd)
  lik_cluster <- apply(lik_cluster,1,function(x,w)return(prod(x^w)),w=w)
  return(lik_cluster)
}

#' @import ucminf
ocmmEst <- function(start, v, x, z, weights, link, gfun, method, rnd=NULL, n_nodes, quad, iclusters){
  len_beta <- ncol(x)
  if (gfun == "glf") {
    if (link == "logit"){
      if (method == "optim"){
        fit <- optim(par=start,negloglik_glf_rnd, v=v, d.matrix=x, rnd.matrix=z, wts=weights, len_beta=len_beta, rnd=rnd, n_nodes=n_nodes, quad=quad, iclusters=iclusters, method="BFGS", hessian = T)
      } else if (method == "ucminf") {
        fit <- ucminf(par=start,negloglik_glf_rnd, v=v, d.matrix=x, rnd.matrix=z, wts=weights, len_beta=len_beta, rnd=rnd, n_nodes=n_nodes, quad=quad, iclusters=iclusters, hessian = 3)
      } else {
        stop("Optimization method not implemented.")
      }
    } else {
      stop("link function not implemented.")
    }
  } else {
    stop("g function not implemented.")
  }
  H=fit$hessian
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
  sigma_rnd <- coef[len_beta+4]
  
  ## degrees of freedom and standard deviation of residuals
  df <- nrow(x)-ncol(x)-length(par_g)
  fitted.values <- inv.logit(g_glf(v, par_g) + x%*%beta)
  sigma2 <- sum((v - fitted.values)^2)/df
  
  ## compute sigma^2 * (x'x)^-1
  colnames(vcov) <- rownames(vcov) <- c(colnames(x),"M", "B", "T", "sigma_rnd")
  list(coefficients = coef,
       vcov = vcov,
##       sigma = sqrt(sigma2),
       sigma_rnd = sigma_rnd,
       df = df,
       logLik = -fit$value)
}
