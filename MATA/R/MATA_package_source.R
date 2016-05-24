


#' Model-Averaged Tail Area Wald (MATA-Wald) Confidence Interval
#'
#' A function for computing the Model-Averaged Tail Area Wald (MATA-Wald)
#' confidence interval, constructed using single-model estimators
#' and model weights.
#'
#' @details
#' \code{mata.wald} may be used to construct model-averaged 
#' confidence intervals, using the Model-Averaged Tail Area (MATA) construction
#' (see Turek and Fletcher (2012) for details).
#' The idea underlying this construction is similar to that of a model-averaged
#' Bayesian credible interval.  This function returns the lower and upper
#' confidence limits of a MATA-Wald interval.
#' 
#' Two usages are supported.  For the normal linear model,
#' using option \code{normal.lm = TRUE} generates a MATA-Wald confidence interval
#' corresponding to the solutions of equations (2) and (3) of Turek and 
#' Fletcher (2012).  The argument \code{residual.dfs} is required for this usage.
#'
#' When the sampling distribution for the estimator is asymptotically 
#' normal (e.g. MLEs), possibly after a transformation, use option \code{normal.lm = FALSE}.
#' This generates a MATA-Wald confidence interval, possibly
#' on a transformed scale, where back-transformation of both confidence limits 
#' may be necessary.  This corresponds to solutions to the equations in Section 3.2 
#' of Turek and Fletcher (2012).
#' 
#' @aliases mata.wald tailarea.z tailarea.t
#' 
#' @param theta.hats A numeric vector containing the parameter estimates under each 
#'                   candidate model.
#'
#' @param se.theta.hats A numeric vector containing the estimated standard error of each 
#'                      value in \code{theta.hats}.
#'
#' @param model.weights A vector containing the model weights for each candidate 
#'                      model.  Calculated from an information criterion,
#'                      such as AIC or BIC. 
#'                      All model weights must be non-negative, and sum to one.
#'
#' @param alpha The desired lower and upper error rate.  The value 0.025
#'              corresponds to a 95\% MATA-Wald confidence interval, and
#'              0.05 to a 90\% interval. Must be between 0 and 0.5.
#'              Default value is 0.025.
#' 
#' @param normal.lm Logical.  TRUE for the normal linear model case, and 
#'                  FALSE otherwise.  When TRUE, the argument 
#'                  \code{residual.dfs} must also be supplied.
#' 
#' @param residual.dfs A vector containing the residual (error) degrees of freedom 
#'                     under each candidate model.  This argument must be provided 
#'                     when \code{normal.lm = TRUE}.
#' 
#' @author Daniel Turek
#'
#' @export
#' 
#' @references Turek, D. and Fletcher, D. (2012). Model-Averaged Wald Confidence Intervals. Computational Statistics and Data Analysis, 56(9), p.2809-2815.
#' 
#' @examples
#'# Normal linear prediction:
#'# Generate single-model Wald and model-averaged MATA-Wald 95% confidence intervals
#'#
#'# Data 'y', covariates 'x1' and 'x2', all vectors of length 'n'.
#'# 'y' taken to have a normal distribution.
#'# 'x1' specifies treatment/group (factor).
#'# 'x2' a continuous covariate.
#'#
#'# Take the quantity of interest (theta) as the predicted response 
#'# (expectation of y) when x1=1 (second group/treatment), and x2=15.
#'
#'n = 20                              # 'n' is assumed to be even
#'x1 = c(rep(0,n/2), rep(1,n/2))      # two groups: x1=0, and x1=1
#'x2 = rnorm(n, mean=10, sd=3)
#'y = rnorm(n, mean = 3*x1 + 0.1*x2)  # data generation
#'
#'x1 = factor(x1)
#'m1 = glm(y ~ x1)                    # using 'glm' provides AIC values.
#'m2 = glm(y ~ x1 + x2)               # using 'lm' doesn't.
#'aic = c(m1$aic, m2$aic)
#'delta.aic = aic - min(aic)
#'model.weights = exp(-0.5*delta.aic) / sum(exp(-0.5*delta.aic))
#'residual.dfs = c(m1$df.residual, m2$df.residual)
#'
#'p1 = predict(m1, se=TRUE, newdata=list(x1=factor(1), x2=15))
#'p2 = predict(m2, se=TRUE, newdata=list(x1=factor(1), x2=15))
#'theta.hats = c(p1$fit, p2$fit)
#'se.theta.hats = c(p1$se.fit, p2$se.fit)
#'
#'#  AIC model weights
#'model.weights
#'
#'#  95% Wald confidence interval for theta (under Model 1)
#'theta.hats[1] + c(-1,1)*qt(0.975, residual.dfs[1])*se.theta.hats[1]
#'
#'#  95% Wald confidence interval for theta (under Model 2)
#'theta.hats[2] + c(-1,1)*qt(0.975, residual.dfs[2])*se.theta.hats[2]
#'
#'#  95% MATA-Wald confidence interval for theta (model-averaging)
#'mata.wald(theta.hats=theta.hats, se.theta.hats=se.theta.hats, 
#'         model.weights=model.weights, normal.lm=TRUE, residual.dfs=residual.dfs)
mata.wald = function(theta.hats, se.theta.hats, model.weights, normal.lm, residual.dfs, alpha=0.025) {
    if(length(theta.hats) != length(se.theta.hats))    stop('dimension mismatch in arguments')
    if(length(theta.hats) != length(model.weights))    stop('dimension mismatch in arguments')
    if(any(se.theta.hats <= 0))                        stop('negative se.theta.hats')
    if(any(model.weights < 0))                         stop('negative model.weights')
    if(abs(sum(model.weights)-1) > 0.001)              stop('model.weights do not sum to 1')
    if(!is.logical(normal.lm))                         stop('normal.lm must be logical (T/F)')
    if((alpha<=0) | (alpha>=0.5))                      stop('alpha outside of meaningful range')
    
    if(normal.lm) {
        if(missing(residual.dfs))                        stop('must specify residual.dfs when normal.lm = TRUE')
        if(length(theta.hats) != length(residual.dfs))   stop('dimension mismatch in arguments')
        if(any(residual.dfs <= 0))                       stop('negative residual.dfs')
        if(any(residual.dfs != round(residual.dfs)))     stop('non-integer residual.dfs')
        
        theta.L = uniroot(f=tailarea.t, interval=c(-1e10, 1e10),
            theta.hats=theta.hats, se.theta.hats=se.theta.hats, 
            model.weights=model.weights, alpha=alpha, 
            residual.dfs=residual.dfs, tol=1e-10)$root
        theta.U = uniroot(f=tailarea.t, interval=c(-1e10, 1e10),
            theta.hats=theta.hats, se.theta.hats=se.theta.hats, 
            model.weights=model.weights, alpha=1-alpha, 
            residual.dfs=residual.dfs, tol=1e-10)$root
    }
    if(!normal.lm) {
        theta.L = uniroot(f=tailarea.z, interval=c(-1e10, 1e10), 
            theta.hats=theta.hats, se.theta.hats=se.theta.hats, 
            model.weights=model.weights, alpha=alpha, tol=1e-10)$root
        theta.U = uniroot(f=tailarea.z, interval=c(-1e10, 1e10), 
            theta.hats=theta.hats, se.theta.hats=se.theta.hats, 
            model.weights=model.weights, alpha=1-alpha, tol=1e-10)$root
    }
    c(theta.L, theta.U)
}


tailarea.z = function(theta, theta.hats, se.theta.hats, model.weights, alpha) {
    z.quantiles = (theta-theta.hats)/se.theta.hats
    tailarea = sum(model.weights*pnorm(z.quantiles)) - alpha
    tailarea
}


tailarea.t = function(theta, theta.hats, se.theta.hats, model.weights, alpha, residual.dfs) {
    t.quantiles = (theta-theta.hats)/se.theta.hats
    tailarea = sum(model.weights*pt(t.quantiles, df=residual.dfs)) - alpha
    tailarea
}





