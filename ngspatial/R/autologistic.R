
#' Return a perfect sample from a centered autologistic model.
#'
#' @details This function implements a perfect sampler for the centered autologistic model. The sampler employs coupling from the past. 
#' @param X the design matrix.
#' @param A the adjacency matrix for the underlying graph.
#' @param theta the vector of parameter values: \eqn{\theta = (\beta^\prime, \eta)^\prime}{\theta = (\beta', \eta)'}.
#' @return A vector that is distributed exactly according to the centered autologistic model with the given design matrix and parameter values.
#' @references
#' Moller, J. (1999) Perfect simulation of conditionally specified models. \emph{Journal of the Royal Statistical Society, Series B, Methodological}, \bold{61}, 251--264.
#' @references
#' Propp, J. G. and Wilson, D. B. (1996) Exact sampling with coupled Markov chains and applications to statistical mechanics. \emph{Random Structures and Algorithms}, \bold{9}(1-2), 223--252.
#' @export
#'  @examples \dontrun{ 
#'
#' # Use the 20 x 20 square lattice as the underlying graph.
#'
#' n = 20
#' A = adjacency.matrix(n)
#'
#' # Assign coordinates to each vertex such that the coordinates are restricted to the unit square
#' # centered at the origin.
#'
#' x = rep(0:(n - 1) / (n - 1), times = n) - 0.5
#' y = rep(0:(n - 1) / (n - 1), each = n) - 0.5
#' X = cbind(x, y)                                 # Use the vertex locations as spatial covariates.
#' beta = c(2, 2)                                  # These are the regression coefficients.
#' mu = exp(X %*% beta)
#' mu = mu / (1 + mu)                              # Compute the independence expectations.
#' col1 = "white"
#' col2 = "black"
#' colfunc = colorRampPalette(c(col1, col2))
#'
#' # Now produce a level plot of the independence expectations. This shows that the large-scale
#' # structure corresponds to a probability gradient that increases as one moves from southwest
#' # to northeast.
#'
#' dev.new()
#' lattice::levelplot(mu ~ x * y, aspect = "iso", col.regions = colfunc(n^2))
#'
#' # Simulate a dataset with the above mentioned regression component and eta equal to 0.6. This
#' # value of eta corresponds to dependence that is moderate in strength.
#'
#' theta = c(beta, eta = 0.6)
#' set.seed(123456)
#' Z = rautologistic(X, A, theta)
#'
#' # Create a level plot of the simulated data.
#'
#' dev.new()
#' lattice::levelplot(Z ~ x * y, aspect = "iso", col.regions = c("white", "black"), colorkey = FALSE)
#' } 

rautologistic = function(X, A, theta)
{
    if (missing(X) || ! is.matrix(X))
        stop("You must supply a design matrix.")
    if (missing(A) || ! is.matrix(A) || ! isSymmetric(A) || ! (A == 0 || A == 1))
        stop("You must supply a symmetric binary adjacency matrix.")
    diag(A) = 0
    if (nrow(X) != nrow(A))
        stop("The supplied design matrix and adjacency matrix are not conformable.")
    if (missing(theta) || ! is.vector(theta, mode = "numeric"))
        stop("You must supply a numeric vector of parameter values.")
    p = length(theta)
    if (ncol(X) != p - 1)
        stop("The given design matrix and vector of regression coefficients are not conformable.")
    as.vector(perfsampler$rautologistic(X, A, theta))
}

autologistic.bmse = function(mat)
{
    if (sum(is.na(mat)) > 0)
    {
        warning("The sample contains NAs.")
        bmse = NA
    }
    else
    {
        mat = as.matrix(mat)
        p = ncol(mat)
        bmse = numeric(p)
        for (i in 1:p)
        {
            temp = perfsampler$bmse(mat[, i])
            if (temp == -1) 
                temp = NA
            bmse[i] = temp
        }
    }
    bmse
}

autologistic.logPL = function(theta, X, A, Z)
{
    p = length(theta)
    beta = theta[-p]
    eta = theta[p]
    Xbeta = X %*% beta
    mu = exp(Xbeta) / (1 + exp(Xbeta))
    logPL = Xbeta + eta * A %*% (Z - mu)
    logPL = t(Z) %*% logPL - sum(log(1 + exp(logPL)))
    -logPL
}

autologistic.control = function(method, control, verbose, p)
{
    if (method == "PL")
    {
        if (is.null(control$confint) || length(control$confint) > 1 || ! control$confint %in% c("bootstrap", "none", "sandwich"))
        {
            if (verbose)
                cat("\nControl parameter 'confint' must be one of \"bootstrap\", \"none\", or \"sandwich\". Setting it to the default value of \"sandwich\".\n")
            control$confint = "sandwich"
        }
        if (control$confint == "none")
            return(list(confint = "none"))
        bootit = control$bootit
        if (is.null(bootit) || ! is.numeric(bootit) || length(bootit) > 1 || ! is.wholenumber(bootit) || bootit < 1)
        {
            if (verbose)
                cat("\nControl parameter 'bootit' must be a positive whole number. Setting it to the default value of 1,000.\n")
            control$bootit = 1000
        }
        if (is.null(control$parallel) || ! is.logical(control$parallel) || length(control$parallel) > 1)
        {
            if (verbose)
                cat("\nControl parameter 'parallel' must be a logical value. Setting it to the default value of TRUE.\n")
            control$parallel = TRUE
        }
        if (control$parallel)
        {
            if (requireNamespace("parallel", quietly = TRUE))
            {
                if (is.null(control$type) || length(control$type) > 1 || ! control$type %in% c("SOCK", "PVM", "MPI", "NWS"))
                {
                    if (verbose)
                        cat("\nControl parameter 'type' must be \"SOCK\", \"PVM\", \"MPI\", or \"NWS\". Setting it to \"SOCK\".\n")
                    control$type = "SOCK"
                }
                nodes = control$nodes
                if (is.null(control$nodes) || ! is.numeric(nodes) || length(nodes) > 1 || ! is.wholenumber(nodes) || nodes < 2)
                    stop("Control parameter 'nodes' must be a whole number greater than 2.")
            }
            else
            {
                if (verbose)
                    cat("\nParallel computation requires package parallel. Setting control parameter 'parallel' to FALSE.\n")
                control$parallel = FALSE    
            }
        }
        old = control
        control = list()
        control$confint = old$confint
        control$bootit = old$bootit
        control$parallel = old$parallel
        if (control$parallel)
        {
            control$type = old$type
            control$nodes = old$nodes
        }
    }
    else
    {
        trainit = control$trainit
        if (is.null(trainit) || ! is.numeric(trainit) || length(trainit) > 1 || ! is.wholenumber(trainit) || trainit < 0)
        {
            if (verbose)
                cat("\nControl parameter 'trainit' must be a positive whole number. Setting it to the default value of 100,000.\n")
            control$trainit = 100000
        }
        if (control$trainit < 10000 && verbose)
             warning("Consider increasing the value of 'trainit'.")
        tol = control$tol
        if (is.null(tol) || ! is.numeric(tol) || length(tol) > 1 || tol <= 0)
        {
            if (verbose)
                cat("\nControl parameter 'tol' must be a positive number. Setting it to the default value of 0.01.\n")
            control$tol = 0.01
        }
        minit = control$minit
        if (is.null(minit) || ! is.numeric(minit) || length(minit) > 1 || ! is.wholenumber(minit) || minit < 1)
        {
            if (verbose)
                cat("\nControl parameter 'minit' must be a positive whole number. Setting it to the default value of 10,000.\n")
            control$minit = 10000
        }
        if (control$minit < 10000 && verbose)
            warning("Consider increasing the value of 'minit'.")
        maxit = control$maxit
        if (is.null(maxit) || ! is.numeric(maxit) || length(maxit) > 1 || ! is.wholenumber(maxit) || maxit < control$minit)
        {
            if (verbose)
                cat("\nControl parameter 'maxit' must be a positive whole number greater than or equal to 'minit'.")
            if (control$minit < 1000000)
                control$maxit = 1000000
            else
                control$maxit = control$minit
            if (verbose)
                cat(" Setting it to ", control$maxit, ".\n", sep = "")
        }
        sigma = control$sigma
        if (is.null(sigma) || ! is.numeric(sigma) || length(sigma) > 1 || sigma <= 0)
        {
            if (verbose)
                cat("\nControl parameter 'sigma' must be a positive number. Setting it to the default value of 1,000.\n")
            control$sigma = 1000
        }
        eta.max = control$eta.max
        if (is.null(eta.max) || ! is.numeric(eta.max) || length(eta.max) > 1 || eta.max <= 0)
        {
            if (verbose)
                cat("\nControl parameter 'eta.max' must be a positive number. Setting it to the default value of 2.\n")
            control$eta.max = 2
        }
        if (control$eta.max < 1 && verbose)
            warning("Consider increasing the value of 'eta.max'.")
        old = control
        control = list()
        control$trainit = old$trainit
        control$tol = old$tol
        control$minit = old$minit
        control$maxit = old$maxit
        control$sigma = old$sigma
        control$eta.max = old$eta.max
    }
    control
}

autologistic.fit = function(X, A, Z, hessian)
{
    start = glm(Z ~ X - 1, family = binomial)$coef
    opt = try(optim(c(start, 1), autologistic.logPL, autologistic.grad, X = X, A = A, Z = Z,
              method = "BFGS", control = list(maxit = 1000), hessian = hessian),
              silent = TRUE)
    if (class(opt) == "try-error")
    {
        coefficients = NULL
        fitted.values = NULL
        linear.predictors = NULL
        residuals = NULL
        convergence = NULL
        message = opt[1]
        value = NULL
        I.hat = NULL
    }
    else
    {
        convergence = opt$convergence
        if (convergence != 0)
            message = opt$message
        else
            message = NULL
        p = ncol(X) + 1
        coefficients = opt$par
        names(coefficients)[p] = "eta"
        Xbeta = X %*% coefficients[-p]
        mu = exp(Xbeta)
        mu = mu / (1 + mu)
        autocovariate = A %*% (Z - mu)
        linear.predictors = Xbeta + coefficients[p] * autocovariate
        fitted.values = exp(linear.predictors)
        fitted.values = fitted.values / (1 + fitted.values)
        residuals = Z - fitted.values
        value = opt$value
        if (hessian)
            I.hat = opt$hessian
    }
    object = list(coefficients = coefficients, fitted.values = fitted.values,
                  linear.predictors = linear.predictors, residuals = residuals, 
                  convergence = convergence, message = message, value = value)
    if (hessian && ! is.null(I.hat))
        object$I.hat = I.hat
    class(object) = "autologistic"
    object
}

autologistic.bootstrap.helper = function(dummy, X, A, theta)
{
    Z = rautologistic(X, A, theta)
    fit = autologistic.fit(X, A, Z, FALSE)
    fit
}

autologistic.sandwich.helper = function(dummy, X, A, theta)
{
    Z = rautologistic(X, A, theta)
    len = length(theta)
    mu = as.vector(exp(X %*% theta[-len]))
    mu = mu / (1 + mu)
    p = as.vector(exp(X %*% theta[-len] + theta[len] * A %*% (Z - mu)))
    p = p / (1 + p)
    c((Z - p) %*% (X - theta[len] * A %*% (X * mu * (1 - mu))), (Z - p) %*% A %*% (Z - mu))
}

autologistic.grad = function(params, X, A, Z)
{
    len = length(params)
    mu = as.vector(exp(X %*% params[-len]))
    mu = mu / (1 + mu)
    p = as.vector(exp(X %*% params[-len] + params[len] * A %*% (Z - mu)))
    p = p / (1 + p)
    -c((Z - p) %*% (X - params[len] * A %*% (X * mu * (1 - mu))), (Z - p) %*% A %*% (Z - mu))
}

autologistic.sandwich = function(X, A, theta, type, bootit, parallel, nodes)
{
    meat = 0
    if (! parallel)
    {
        for (j in 1:bootit)
        {
            Z = rautologistic(X, A, theta)
            len = length(theta)
            mu = as.vector(exp(X %*% theta[-len]))
            mu = mu / (1 + mu)
            p = as.vector(exp(X %*% theta[-len] + theta[len] * A %*% (Z - mu)))
            p = p / (1 + p)
            gr = c((Z - p) %*% (X - theta[len] * A %*% (X * mu * (1 - mu))), (Z - p) %*% A %*% (Z - mu))
            meat = meat + gr %o% gr / bootit
        }
    }
    else
    {
		if (requireNamespace("parallel", quietly = TRUE))
		{
            cl = parallel::makeCluster(nodes, type)
            parallel::clusterSetRNGStream(cl, NULL)
            parallel::clusterEvalQ(cl, library(ngspatial))
            gathered = parallel::clusterApplyLB(cl, 1:bootit, autologistic.sandwich.helper, X, A, theta)
            parallel::stopCluster(cl)
            for (j in 1:bootit)
            {
                gr = gathered[[j]]
                meat = meat + gr %o% gr / bootit                
            }
		}
    }
    meat
}

autologistic.bootstrap = function(X, A, theta, type, bootit, parallel, nodes)
{
    boot.sample = data.frame(matrix(, bootit, length(theta)))
    if (! parallel)
    {
        for (j in 1:bootit)
        {
            fit = autologistic.bootstrap.helper(NULL, X, A, theta)
            temp = rep(NA, length(theta))
            if (is.null(fit$convergence) || fit$convergence != 0)
                warning(fit$message)
            else
                temp = fit$coef
            boot.sample[j, ] = temp
        }
    }
    else
    {
		if (requireNamespace("parallel", quietly = TRUE))
		{
            cl = parallel::makeCluster(nodes, type)
            parallel::clusterSetRNGStream(cl, NULL)
            parallel::clusterEvalQ(cl, library(ngspatial))
            gathered = parallel::clusterApplyLB(cl, 1:bootit, autologistic.bootstrap.helper, X, A, theta)
            parallel::stopCluster(cl)
            for (j in 1:bootit)
            {
                fit = gathered[[j]]
                temp = rep(NA, length(theta))
                if (is.null(fit$convergence) || fit$convergence != 0)
                    warning(fit$message)
                else
                    temp = fit$coef
                boot.sample[j, ] = temp
			}         
        }
    }
    boot.sample
}

#' Extract model residuals.
#'
#' @param object an object of class \code{autologistic}, typically the result of a call to \code{\link{autologistic}}.
#' @param type the type of residuals that should be returned. The alternatives are \dQuote{\code{deviance}} (default), \dQuote{\code{pearson}}, and \dQuote{\code{response}}.
#' @param \dots additional arguments.
#' @return A vector of residuals.
#' @seealso \code{\link{autologistic}}, \code{\link{residuals.glm}}
#' @method residuals autologistic
#' @export

residuals.autologistic = function(object, type = c("deviance", "pearson", "response"), ...)
{
    type = match.arg(type)
    if (type == "response")
        return(object$residuals)
    else if (type == "deviance")
    {
        if (is.null(object$y))
            y = object$residuals + object$fitted.values
        phat = object$fitted.values
        d = numeric(length(y))
        zero = which(y == 0)
        d[zero] = -2 * log(1 - phat[zero])
        one = which(y == 1)
        d[one] = -2 * log(phat[one])
        return(sqrt(d) * sign(object$residuals))
    }
    else # type == "pearson"
    {
        phat = object$fitted.values
        se = sqrt(phat * (1 - phat))
        return(object$residuals / se)
    }
}

#' Return the estimated covariance matrix for an \code{autologistic} model object.
#'
#' @param object a fitted \code{autologistic} model object.
#' @param \dots additional arguments.
#' @return An estimate of the covariance matrix of the parameters (in a Bayesian setting), or an estimate of the covariance matrix of the maximum pseudolikelihood estimator of the parameters.
#' @method vcov autologistic
#' @export

vcov.autologistic = function(object, ...)
{
    if (is.null(object$sample))
        V = object$S
    else if (sum(is.na(object$sample)) > 0)
        stop("The sample contains NAs.")
    else
    {
        V = cov(object$sample)
        rownames(V) = colnames(V) = names(object$coef)
    }
    V
}

#' Print a summary of a centered autologistic model fit.
#'
#' @details This function displays (1) the call to \code{\link{autologistic}}, (2) the values of the control parameters, (3) a table of estimates, and (4) the size of the bootstrap/posterior sample. Each row of the table of estimates shows a parameter estimate, the confidence or HPD interval for the parameter, and, where applicable, the Monte Carlo standard error.
#' @param object an object of class \code{autologistic}, typically the result of a call to \code{\link{autologistic}}.
#' @param alpha the significance level for the quantile/HPD intervals. The default is 0.05.
#' @param digits the number of significant digits to display. The default is 4.
#' @param \dots additional arguments.
#' @seealso \code{\link{autologistic}}
#' @method summary autologistic
#' @export

summary.autologistic = function(object, alpha = 0.05, digits = 4, ...)
{
	if (! is.numeric(alpha) || length(alpha) > 1 || alpha <= 0 || alpha >= 1)
	    stop("'alpha' must be a number between 0 and 1.")
    cat("\nCall:\n\n")
    print(object$call)
    cat("\nControl parameters:\n")
    control.table = cbind(unlist(object$control))
    colnames(control.table) = ""
    print(control.table, quote = FALSE)
    p = length(object$coef)
    if (is.null(object$sample))
    {
        if (is.null(object$S))
            coef.table = cbind(object$coef, NA, NA, NA)
        else
        {
            scale = qnorm(1 - alpha / 2)
            se = sqrt(diag(object$S))
            ci = cbind(object$coef - scale * se, object$coef + scale * se)
            coef.table = cbind(object$coef, ci, NA)
        }
    }
    else if (sum(is.na(object$sample)) > 0)
    {
        warning("The sample contains NAs.")
        coef.table = cbind(object$coef, NA, NA, NA)
    }
    else
    {
        ci = matrix(, p, 2)
        if (object$method == "PL")
        {
            for (j in 1:p)
                ci[j, ] = quantile(object$sample[, j], c(alpha / 2, 1 - alpha / 2), na.rm = TRUE)
        }
        else if (object$method == "Bayes")
        {
            for (j in 1:p)
                ci[j, ] = hpd(object$sample[, j], alpha)
        }
        coef.table = cbind(object$coef, ci, object$mcse)
    }
    colnames(coef.table) = c("Estimate", "Lower", "Upper", "MCSE")
    rownames(coef.table) = names(object$coef)
    cat("\nCoefficients:\n\n")
    print(signif(coef.table, digits))
    cat("\nNumber of iterations:", object$iter, "\n\n")
}

#' Fit a centered autologistic model using maximum pseudolikelihood estimation or MCMC for Bayesian inference.
#'
#' @details This function fits the centered autologistic model of Caragea and Kaiser (2009) using maximum pseudolikelihood estimation or MCMC for Bayesian inference. 
#'          The joint distribution for the centered autologistic model is 
#'          \deqn{\pi(Z\mid\theta)=c(\theta)^{-1}\exp\left(Z^\prime X\beta - \eta Z^\prime A\mu + \frac{\eta}{2}Z^\prime AZ\right),}{\pi(Z | \theta)=c(\theta)^{-1} exp(Z'X\beta - \eta Z'A\mu + 0.5 \eta Z'AZ),}
#'          where \eqn{\theta = (\beta^\prime, \eta)^\prime}{\theta = (\beta', \eta)'} is the parameter vector, \eqn{c(\theta)} is an intractable normalizing function, \eqn{Z} is the response vector, \eqn{X} is the design matrix, 
#'          \eqn{\beta} is a \eqn{(p-1)}-vector of regression coefficients, \eqn{A} is the adjacency matrix for the underlying graph, \eqn{\mu} is the vector of independence expectations, 
#'          and \eqn{\eta} is the spatial dependence parameter. 
#'          \cr
#'          \cr
#'          Maximum pseudolikelihood estimation sidesteps the intractability of \eqn{c(\theta)} by maximizing the product of the conditional likelihoods.
#'          Confidence intervals are then obtained using a parametric bootstrap or a normal approximation, i.e., sandwich estimation. The bootstrap datasets are generated by perfect sampling (\code{\link{rautologistic}}).
#'          The bootstrap samples can be generated in parallel using the \pkg{snow} package.
#'          \cr
#'          \cr
#'          Bayesian inference is obtained using the auxiliary variable algorithm of Moller et al. (2006).
#'          The auxiliary variables are generated by perfect sampling.
#'          \cr
#'          \cr
#'          The prior distributions are (1) zero-mean normal with independent coordinates for \eqn{\beta}, and (2) uniform for \eqn{\eta}.
#'          The common standard deviation for the normal prior can be supplied by the user as control parameter \code{sigma}. The default is 1,000. The uniform prior has support [0, 2] by default, but the right endpoint can be supplied (as control parameter \code{eta.max}) by the user. 
#'          \cr
#'          \cr
#'          The posterior covariance matrix of \eqn{\theta} is estimated using samples obtained during a training run. The default number of iterations for the training run is 100,000, but this can be controlled by the user (via parameter \code{trainit}). The estimated covariance matrix is then used as the proposal variance for a Metropolis-Hastings random walk. The proposal distribution is normal. The posterior samples obtained during the second run are used for inference. The length of the run can be controlled by the user via parameters \code{minit}, \code{maxit}, and \code{tol}. The first determines the minimum number of iterations. If \code{minit} has been reached, the sampler will terminate when \code{maxit} is reached or all Monte Carlo standard errors are smaller than \code{tol}, whichever happens first. The default values for \code{minit}, \code{maxit}, and \code{tol} are 100,000; 1,000,000; and 0.01, respectively.
#'
#' @param formula an object of class \code{\link{formula}}: a symbolic description of the model to be fitted.
#' @param data an optional data frame, list, or environment (or object coercible by \code{\link{as.data.frame}} to a data frame) containing the variables in the model. If not found in \code{data}, the variables are taken from \code{environment(formula)}, typically the environment from which \code{autologistic} is called.
#' @param A the adjacency matrix for the underlying graph.
#' @param method the method to use for inference. \dQuote{\code{PL}} (the default) enables maximum pseudolikelihood estimation, and \dQuote{\code{Bayes}} enables Bayesian inference.
#' @param control a list of the following control parameters.
#'        \describe{
#'        \item{\code{confint}}{ (PL) the method for computing confidence intervals. The possible values are \dQuote{\code{sandwich}} (the default), \dQuote{\code{bootstrap}}, and \dQuote{\code{none}}.}
#'        \item{\code{sigma}}{ (Bayes) the common standard deviation for \eqn{\beta}'s prior distribution. Defaults to 1,000.}
#'        \item{\code{eta.max}}{ (Bayes) the right endpoint for \eqn{\eta}'s prior distribution. Defaults to 2.}
#'        \item{\code{bootit}}{ (PL) the size of the bootstrap sample. This applies when \code{confint} is \dQuote{\code{sandwich}} or \dQuote{\code{bootstrap}}, since samples from the fitted model are needed in both cases. Defaults to 1,000.}
#'        \item{\code{trainit}}{ (Bayes) the length of the training run. Defaults to 100,000.}
#'        \item{\code{minit}}{ (Bayes) the minimum number of posterior samples. Defaults to 100,000.}
#'        \item{\code{maxit}}{ (Bayes) the maximum number of posterior samples. Defaults to 1,000,000.}
#'        \item{\code{tol}}{ (Bayes) the tolerance. Defaults to 0.01.}
#'        \item{\code{parallel}}{ (PL) a logical value indicating whether to parallelize the bootstrap. This defaults to \code{TRUE} if the \pkg{snow} package is installed.}
#'        \item{\code{type}}{ (PL) the cluster type, one of \dQuote{\code{SOCK}} (default), \dQuote{\code{PVM}}, \dQuote{\code{MPI}}, or \dQuote{\code{NWS}}.}
#'        \item{\code{nodes}}{ (PL) the number of slave nodes to create.}}
#' @param model a logical value indicating whether the model frame should be included as a component of the returned value.
#' @param x a logical value indicating whether the model matrix used in the fitting process should be returned as a component of the returned value.
#' @param y a logical value indicating whether the response vector used in the fitting process should be returned as a component of the returned value.
#' @param verbose a logical value indicating whether to print various messages to the screen, including progress updates when \code{method} is \dQuote{\code{Bayes}}. Defaults to \code{FALSE}.
#' @return \code{autologistic} returns an object of class \dQuote{\code{autologistic}}, which is a list containing the following components.
#'         \item{coefficients}{the point estimate of \eqn{\theta}.}
#'         \item{fitted.values}{the fitted mean values, obtained by transforming the linear predictors by the inverse of the link function.}
#'         \item{linear.predictors}{the linear fit on link scale.}
#'         \item{residuals}{the response residuals.}
#'         \item{iter}{the size of the bootstrap/posterior sample.}
#'         \item{sample}{(where relevant) an \code{iter} by \eqn{p} matrix containing the bootstrap/posterior samples.}
#'         \item{mcse}{(where relevant) a \eqn{p}-vector of Monte Carlo standard errors.}
#'         \item{S}{(where relevant) the estimated sandwich matrix.}
#'         \item{accept}{(Bayes) the acceptance rate for the MCMC sampler.}
#'         \item{y}{if requested (the default), the \code{y} vector used.}
#'         \item{X}{if requested, the model matrix.}
#'         \item{model}{if requested (the default), the model frame.}
#'         \item{call}{the matched call.}
#'         \item{formula}{the formula supplied.}
#'         \item{method}{the method used for inference.}
#'         \item{convergence}{the integer code returned by \code{\link{optim}} subsequent to optimizing the pseudolikelihood.}
#'         \item{message}{a character string to go along with \code{convergence}.}
#'         \item{terms}{the \code{\link{terms}} object used.}
#'         \item{data}{the \code{data} argument.}
#'         \item{xlevels}{(where relevant) a record of the levels of the factors used in fitting.}
#'         \item{control}{a list containing the names and values of the control parameters.}
#'         \item{value}{the value of the negative log pseudolikelihood at its minimum.}
#' @references
#' Caragea, P. and Kaiser, M. (2009) Autologistic models with interpretable parameters. \emph{Journal of Agricultural, Biological, and Environmental Statistics}, \bold{14}(3), 281--300.
#' @references
#' Hughes, J., Haran, M. and Caragea, P. C. (2011) Autologistic models for binary data on a lattice. \emph{Environmetrics}, \bold{22}(7), 857--871. 
#' @references
#' Moller, J., Pettitt, A., Berthelsen, K., and Reeves, R. (2006) An efficient Markov chain Monte Carlo method for distributions with intractable normalising constants. \emph{Biometrika}, \bold{93}(2), 451--458.
#' @seealso \code{\link{rautologistic}}, \code{\link{residuals.autologistic}}, \code{\link{summary.autologistic}}, \code{\link{vcov.autologistic}}
#' @export
#'  @examples \dontrun{
#'
#' # Use the 20 x 20 square lattice as the underlying graph.
#'
#' n = 20
#' A = adjacency.matrix(n)
#'
#' # Assign coordinates to each vertex such that the coordinates are restricted to the unit square
#' # centered at the origin.
#'
#' x = rep(0:(n - 1) / (n - 1), times = n) - 0.5
#' y = rep(0:(n - 1) / (n - 1), each = n) - 0.5
#' X = cbind(x, y)                                 # Use the vertex locations as spatial covariates.
#' beta = c(2, 2)                                  # These are the regression coefficients.
#' col1 = "white"
#' col2 = "black"
#' colfunc = colorRampPalette(c(col1, col2))
#'
#' # Simulate a dataset with the above mentioned regression component and eta equal to 0.6. This
#' # value of eta corresponds to dependence that is moderate in strength.
#'
#' theta = c(beta, eta = 0.6)
#' set.seed(123456)
#' Z = rautologistic(X, A, theta)
#'
#' # Find the MPLE, and do not compute confidence intervals.
#'
#' fit = autologistic(Z ~ X - 1, A = A, control = list(confint = "none"))
#' summary(fit)
#'
#' # Compute confidence intervals based on the normal approximation. Estimate the "filling" in the
#' # sandwich matrix using a parallel parametric bootstrap, where the computation is distributed
#' # across six cores on the host workstation.
#'
#' set.seed(123456)
#' fit = autologistic(Z ~ X - 1, A = A, verbose = TRUE,
#'                    control = list(confint = "sandwich", nodes = 6))
#' summary(fit)
#'
#' # Compute confidence intervals based on a parallel parametric bootstrap. Use a bootstrap sample
#' # of size 500, and distribute the computation across six cores on the host workstation.
#'
#' set.seed(123456)
#' fit = autologistic(Z ~ X - 1, A = A, verbose = TRUE,
#'                    control = list(confint = "bootstrap", bootit = 500, nodes = 6))
#' summary(fit)
#'
#' # Make some level plots of the residuals.
#'
#' dev.new()
#' lattice::levelplot(residuals(fit) ~ x * y, aspect = "iso", col.regions = colfunc(n^2))
#' dev.new()
#' lattice::levelplot(residuals(fit, type = "pearson") ~ x * y, aspect = "iso",
#'                    col.regions = colfunc(n^2))
#' dev.new()
#' lattice::levelplot(residuals(fit, type = "response") ~ x * y, aspect = "iso",
#'                    col.regions = colfunc(n^2))
#'
#' # Do MCMC for Bayesian inference. The length of the training run will be 10,000, and
#' # the length of the subsequent inferential run will be at least 10,000.
#'
#' set.seed(123456)
#' fit = autologistic(Z ~ X - 1, A = A, verbose = TRUE, method = "Bayes",
#'                    control = list(trainit = 10000, minit = 10000, sigma = 1000))
#' summary(fit)
#' } 

autologistic = function(formula, data, A, method = c("PL", "Bayes"), model = TRUE,
                        x = FALSE, y = FALSE, verbose = FALSE, control = list())
{
    cl = match.call()   
    if (missing(formula))
        stop("You must supply a formula.")
    if (missing(A) || ! is.matrix(A) || ! isSymmetric(A) || ! (A == 0 || A == 1))
        stop("You must supply a symmetric binary adjacency matrix.")
    diag(A) = 0
    mf = match.call(expand.dots = FALSE)
    m = match(c("formula", "data"), names(mf), 0)
    mf = mf[c(1, m)]
    mf[[1]] = as.name("model.frame")
    mf = eval(mf, parent.frame())
    mt = attr(mf, "terms")
    Z = model.response(mf, "numeric")
    X = model.matrix(mt, mf)
    if (! is.vector(Z) || ! (Z == 0 || Z == 1))
        stop("The response must be a binary vector.")
    if (sum(c(nrow(X), nrow(A)) != length(Z)) > 0)
        stop("The supplied response vector/design matrix/adjacency matrix are not conformable.")
    method = match.arg(method)
    if (! is.logical(verbose) || length(verbose) > 1)
        stop("'verbose' must be a logical value.")
    p = ncol(X)
    if (! is.list(control))
        stop("'control' must be a list.")
    control = autologistic.control(method, control, verbose, p)
    hessian = FALSE
    if (method == "PL" && control$confint == "sandwich")
        hessian = TRUE
    fit = autologistic.fit(X, A, Z, hessian)
    if (is.null(fit$convergence) || fit$convergence != 0)
        stop(fit$message, " Maximum pseudolikelihood estimation failed.")
    if (method == "PL")
    {
        fit$iter = 0
        if (! is.null(control$bootit))
            fit$iter = control$bootit
        if (control$confint == "sandwich")
        {
            if (verbose)
            {
                cat("\nWarning: Bootstrapping may be time consuming.\n\n")
                flush.console()
                Sys.sleep(0.5)
            }
            meat = autologistic.sandwich(X, A, fit$coef, control$type, control$bootit, control$parallel, control$nodes)
            bread = try(solve(fit$I.hat), silent = TRUE)
            if (class(bread) != "try-error")
                fit$S = bread %*% meat %*% bread
            else
                fit$S = NULL
            fit = fit[-which(names(fit) == "I.hat")]
            class(fit) = "autologistic"
        }
        else if (control$confint == "bootstrap")
        {
            if (verbose)
            {
                cat("\nWarning: Bootstrapping may be time consuming.\n\n")
                flush.console()
                Sys.sleep(0.5)
            }
            fit$sample = autologistic.bootstrap(X, A, fit$coef, control$type, control$bootit, control$parallel, control$nodes)
            fit$mcse = autologistic.bmse(fit$sample)
        }
    }
    else # method == "Bayes"
    {
        if (verbose)
        {
            cat("\nWarning: MCMC may be time consuming.\n\n")
            flush.console()
            Sys.sleep(0.5)
        }
        s2 = rep(control$sigma^2, p)
        temp = Moller.run(X, A, Z, fit$coef, control$trainit, control$tol, control$minit, control$maxit, s2, control$eta.max, verbose)
        fit = c(fit[-which(names(fit) == "coefficients")], temp)
        class(fit) = c("autologistic")
        p = p + 1
        beta = as.matrix(fit$sample[, -p])
        eta = as.vector(fit$sample[, p])
        n = length(Z)
        linear.predictors = numeric(n)
        fitted.values = numeric(n)
        iter = fit$iter
        for (j in 1:iter)
        {
            Xbeta = X %*% beta[j, ]
            mu = exp(Xbeta)
            mu = mu / (1 + mu)
            auto = A %*% (Z - mu)
            linpred = Xbeta + eta[j] * auto
            linear.predictors = linear.predictors + linpred / iter
            fv = exp(linpred)
            fv = fv / (1 + fv)
            fitted.values = fitted.values + fv / iter
        }
        residuals = Z - fitted.values
        fit$linear.predictors = linear.predictors
        fit$fitted.values = fitted.values
        fit$residuals = residuals
        fit$accept = sum(diff(fit$sample[, 1]) != 0) / fit$iter
    }
    names(fit$coefficients) = c(colnames(X), "eta")
    fit$xlevels = .getXlevels(mt, mf)
    fit$call = cl
    fit$terms = mt
    fit$method = method
    fit$control = control
    if (model)
        fit$model = mf
    if (x)
        fit$x = X
    if (y)
        fit$y = Z
    fit
}

