#' Bayesian inference for multiple linear regression
#' 
#' bayes.lm is used to fit linear models in the Bayesian paradigm. It can be used to carry out regression, 
#' single stratum analysis of variance and analysis of covariance (although these are not tested). This
#' documentation is shamelessly adapated from the lm documentation
#' 
#' 
#' @param formula an object of class \code{\link[stats]{formula}} (or one that can be coerced to that class): a symbolic 
#' description of the model to be fitted. The details of model specification are given under `Details'.
#' @param data an optional data frame, list or environment (or object coercible by \code{\link[base]{as.data.frame}} to a 
#' data frame) containing the variables in the model. If not found in data, the variables are taken 
#' from \code{environment(formula)}, typically the environment from which \code{bayes.lm} is called.
#' @param subset an optional vector specifying a subset of observations to be used in the fitting process.
#' @param na.action a function which indicates what should happen when the data contain \code{NA}s. The 
#' default is set by the \code{\link[stats]{na.action}} setting of options, and is \code{link[stats]{na.fail}}
#'  if that is unset. The `factory-fresh' default is \code{\link[stats]{na.omit}}. Another possible value 
#'  is \code{NULL}, no action. Value \code{\link[stats]{na.exclude}} can be useful.
#' @param model,x,y logicals. If \code{TRUE} the corresponding components of the fit (the model frame, the model matrix, the response)
#'  are returned.
#' \eqn{\beta}{beta}. This argument is ignored for a flat prior.
#' @param center logical or numeric. If \code{TRUE} then the covariates will be centered on their means to make them
#' orthogonal to the intercept. This probably makes no sense for models with factors, and if the argument
#' is numeric then it contains a vector of covariate indices to be centered (not implemented yet).
#' @param prior A list containing b0 (A vector of prior coefficients) and V0 (A prior covariance matrix)
#' @param sigma the population standard deviation of the errors. If \code{FALSE} then this is estimated from the residual sum of squares from the ML fit.
#' 
#' @details Models for \code{bayes.lm} are specified symbolically. A typical model has the form 
#' \code{response ~ terms} where \code{response} is the (numeric) response vector and \code{terms} is a
#'  series of terms which specifies a linear predictor for \code{response}. A terms specification of the 
#'  form \code{first + second} indicates all the terms in \code{first} together with all the terms in 
#'  \code{second} with duplicates removed. A specification of the form \code{first:second} indicates the 
#'  set of terms obtained by taking the interactions of all terms in \code{first} with all terms in 
#'  \code{second}. The specification \code{first*second} indicates the cross of \code{first} and \code{second}.
#'   This is the same as \code{first + second + first:second}.
#'   
#'   See \code{\link[stats]{model.matrix}} for some further details. The terms in the formula will be 
#'   re-ordered so that main effects come first, followed by the interactions, all second-order, 
#'   all third-order and so on: to avoid this pass a \code{terms} object as the formula 
#'   (see \code{\link[stats]{aov}} and \code{demo(glm.vr)} for an example).
#'   
#'   A formula has an implied intercept term. To remove this use either \code{y ~ x - 1} or 
#'   \code{y ~ 0 + x}. See \code{\link[stats]{formula}} for more details of allowed formulae.
#'   
#'   \code{bayes.lm} calls the lower level function \code{lm.fit} to get the maximum likelihood estimates
#'    see below, for the actual numerical computations. For programming only, you may consider doing 
#'    likewise.
#'    
#'    \code{subset} is evaluated in the same way as variables in formula, that is first in data and 
#'    then in the environment of formula.
#'    
#' @return \code{bayes.lm} returns an object of class \code{Bolstad}.
#' The \code{summary} function is used to obtain and print a summary of the results much like the usual 
#' summary from a linear regression using \code{\link[stats]{lm}}.
#' The generic accessor functions \code{coef, fitted.values and residuals}
#' extract various useful features of the value returned by \code{bayes.lm}. Note that the residuals
#' are computed at the posterior mean values of the coefficients.
#' 
#' An object of class "Bolstad" from this function is a list containing at least the following components:
#' \item{coefficients}{a named vector of coefficients which contains the posterior mean}
#' \item{post.var}{a matrix containing the posterior variance-covariance matrix of the coefficients}
#' \item{post.sd}{sigma}
#' \item{residuals}{the residuals, that is response minus fitted values (computed at the posterior mean)}
#' \item{fitted.values}{the fitted mean values (computed at the posterior mean)}
#' \item{df.residual}{the residual degrees of freedom}
#' \item{call}{the matched call}
#' \item{terms}{the \code{\link[stats]{terms}} object used}
#' \item{y}{if requested, the response used}
#' \item{x}{if requested, the model matrix used}
#' \item{model}{if requested (the default), the model frame used}
#' \item{na.action}{(where relevant) information returned by \code{model.frame} on the special 
#' handling of \code{NA}s}
#' 
#' @keywords misc
#' @examples
#' data(bears)
#' bears = subset(bears, Obs.No==1)
#' bears = bears[,-c(1,2,3,11,12)]
#' bears = bears[ ,c(7, 1:6)]
#' bears$Sex = bears$Sex - 1
#' log.bears = data.frame(log.Weight = log(bears$Weight), bears[,2:7])
#' 
#' b0 = rep(0, 7)
#' V0 = diag(rep(1e6,7))
#' 
#' fit = bayes.lm(log(Weight)~Sex+Head.L+Head.W+Neck.G+Length+Chest.G, data = bears,
#'                prior = list(b0 = b0, V0 = V0))
#' summary(fit)
#' print(fit)
#' 
#' 
#' ## Dobson (1990) Page 9: Plant Weight Data:
#' ctl <- c(4.17,5.58,5.18,6.11,4.50,4.61,5.17,4.53,5.33,5.14)
#' trt <- c(4.81,4.17,4.41,3.59,5.87,3.83,6.03,4.89,4.32,4.69)
#' group <- gl(2, 10, 20, labels = c("Ctl","Trt"))
#' weight <- c(ctl, trt)
#' 
#' lm.D9 <- lm(weight ~ group)
#' bayes.D9 <- bayes.lm(weight ~ group)
#' 
#' summary(lm.D9)
#' summary(bayes.D9)
#' 
#' @export bayes.lm

bayes.lm = function(formula, data, subset, na.action, model = TRUE, x = FALSE, y = FALSE,
                    center = TRUE, prior = NULL, sigma = FALSE){
  ret.x = x
  ret.y = y
  cl = match.call()
  mf = match.call(expand.dots = FALSE)
  m = match(c("formula", "data", "subset", "na.action"), names(mf), 0L)
  mf = mf[c(1L, m)]
  mf$drop.unused.levels = TRUE
  mf[[1L]] = quote(stats::model.frame)
  mf = eval(mf, parent.frame())
  mt = attr(mf, "terms")
  y = model.response(mf, "numeric")
  
  if (is.empty.model(mt)) {
#     x = NULL
#     z = list(coefficients = if (is.matrix(y)) matrix(, 0, 
#                                                       3) else numeric(), residuals = y, fitted.values = 0 * 
#                 y, rank = 0L, df.residual = if (!is.null(w)) sum(w != 
#                                                                                 0) else if (is.matrix(y)) nrow(y) else length(y))
  }
  else {
    x = model.matrix(mt, mf, contrasts)
    if(center){
      if(is.logical(center)){
        np = ncol(x)
        x[,2:np] = scale(x[,2:np], scale = FALSE, center = TRUE)
      }
    }
      
    z = z.ls = lm.fit(x, y)
    p1 = 1:z$rank
    z$cov.unscaled = chol2inv(z$qr$qr[p1, p1, drop = FALSE])
    
    z$prior = prior
    if(!is.null(prior)){
      prior.prec = solve(prior$V0)
      resVar = if(is.logical(sigma) && !sigma){
        sum(z$residuals^2) / z$df.residual
      }else{
        sigma^2
      }
      
      ls.prec = solve(resVar * z$cov.unscaled)
      post.prec = prior.prec + ls.prec
      V1 = solve(post.prec)
      b1 = V1 %*% prior.prec %*% prior$b0 + V1 %*% ls.prec %*% coef(z.ls)
      z$post.mean = z$coefficients = as.vector(b1)
      z$post.var = V1
      z$post.sd = sqrt(resVar)
      
    }else{
      resVar = if(is.logical(sigma) && !sigma){
        sum(z$residuals^2) / z$df.residual
      }else{
        sigma^2
      }
      z$post.mean = z$coefficients
      z$post.var = resVar * z$cov.unscaled
      z$post.sd = sqrt(resVar)
    }
    
    z$fitted.values = x %*% z$post.mean
    z$residuals = y - z$fitted.values
    z$df.residual = nrow(x) - ncol(x)
  }
  class(z) = c("Bolstad", "lm")
  z$na.action = attr(mf, "na.action")
  z$call = cl
  z$terms = mt
  if (model) 
    z$model = mf
  if (ret.x) 
    z$x = x
  if (ret.y) 
    z$y = y
  
  
  
  z
}