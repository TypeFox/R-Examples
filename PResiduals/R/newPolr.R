#' slightly modified version of polr from MASS
#' @param formula a formula
#' @param data an optional data frame, list or environment (or object coercible
#' by \code{\link{as.data.frame}} to a data frame) containing the variables in
#' the model.  If not found in \code{data}, the variables are taken from
#' \code{environment(formula)}, typically the environment from which
#' \code{cobot} is called.
#' @param subset an optional vector specifying a subset of observations to be
#' used in the fitting process.
#' @param weights optional case weights in fitting.  Default to 1.
#' @param na.action a function which indicates what should happen when the data
#' contain \code{NA}s.  The default is is \code{\link{na.fail}}.  Another
#' possible value is \code{NULL}, no action.  Value \code{\link{na.exclude}} can
#' be useful.
#' @param offset this can be used to specify an \emph{a priori} known component
#' to be included in the linear predictor during fitting.  This should be
#' \code{NULL} or a numeric vector of length equal to the number of cases.  One
#' or more \code{\link{offset}} terms can be included in the formula instead or
#' as well, and if more than one are specified their sum is used.  See
#' \code{\link{model.offset}}.
#' @param contrasts     a list of contrasts to be used for some or all of the
#' factors appearing as variables in the model formula.
#' @param ...     additional arguments to be passed to \code{\link{optim}}, most
#' often a \code{control} argument.
#' @param Hess logical for whether the Hessian (the observed information matrix)
#' should be returned.  Use this if you intend to call \code{\link{summary}} or
#' \code{\link{vcov}} on the fit.
#' @param model logical for whether the model matrix should be returned.
#' @param method logistic or probit or complementary log-log or cauchit
#' (corresponding to a Cauchy latent variable).
#' @return   A object of class \code{"polr"}.  This has components
#' \item{coefficients}{the coefficients of the linear predictor, which has no
#' intercept.}
#' \item{zeta}{the intercepts for the class boundaries.}
#' \item{deviance}{the residual deviance.}
#' \item{fitted.values}{a matrix, with a column for each level of the response.}
#' \item{lev}{the names of the response levels.}
#' \item{terms}{the \code{terms} structure describing the model.}
#' \item{df.residual}{the number of residual degrees of freedoms, calculated
#' using the weights.}
#' \item{edf}{the (effective) number of degrees of freedom used by the model}
#' \item{n, nobs}{the (effective) number of observations, calculated using the
#' weights. (\code{nobs} is for use by \code{\link{stepAIC}}).}
#' \item{call}{the matched call.}
#' \item{method}{the matched method used.}
#' \item{convergence}{the convergence code returned by \code{optim}.}
#' \item{niter}{the number of function and gradient evaluations used by
#' \code{\link{optim}}.}
#' \item{lp}{the linear predictor (including any offset).}
#' \item{Hessian}{(if \code{Hess} is true).  Note that this is a numerical
#' approximation derived from the optimization proces.}
#' \item{model}{(if \code{model} is true).}
#' @references
#' polr from MASS
#' @seealso \code{\link{optim}}, \code{\link{glm}}, \code{\link[nnet]{multinom}}
#' @importFrom stats model.weights model.offset model.response glm.fit .getXlevels binomial


newpolr <- function(formula, data, weights, start, ..., subset,
                 na.action, contrasts = NULL, Hess = FALSE, model = TRUE,
                 method = c("logit", "probit", "cloglog", "cauchit"))
{
    m <- match.call(expand.dots = FALSE)
    method <- match.arg(method)
    if(is.matrix(eval.parent(m$data))) m$data <- as.data.frame(data)
    m$start <- m$Hess <- m$method <- m$model <- m$... <- NULL
    m[[1L]] <- as.name("model.frame")
    m <- eval.parent(m)
    Terms <- attr(m, "terms")
    x <- model.matrix(Terms, m, contrasts)
    xint <- match("(Intercept)", colnames(x), nomatch = 0L)
    n <- nrow(x)
    pc <- ncol(x)
    cons <- attr(x, "contrasts") # will get dropped by subsetting
    if(xint > 0L) {
        x <- x[, -xint, drop = FALSE]
        pc <- pc - 1L
    } else warning("an intercept is needed and assumed")
    wt <- model.weights(m)
    if(!length(wt)) wt <- rep(1, n)
    offset <- model.offset(m)
    if(length(offset) <= 1L) offset <- rep(0, n)
    y <- model.response(m)
    if(!is.factor(y)) stop("response must be a factor")
    lev <- levels(y); llev <- length(lev)
    if(llev <= 2L) stop("response must have 3 or more levels")
    y <- unclass(y)
    q <- llev - 1L
    Y <- matrix(0, n, q)
    if(missing(start)) {
        # try logistic/probit regression on 'middle' cut
        q1 <- llev %/% 2L
        y1 <- (y > q1)
        X <- cbind(Intercept = rep(1, n), x)
        fit <-
            switch(method,
                   "logit"= glm.fit(X, y1, wt, family = binomial(), offset = offset),
                   "probit" = glm.fit(X, y1, wt, family = binomial("probit"), offset = offset),
                   ## this is deliberate, a better starting point
                   "cloglog" = glm.fit(X, y1, wt, family = binomial("probit"), offset = offset),
                   "cauchit" = glm.fit(X, y1, wt, family = binomial("cauchit"), offset = offset))
        if(!fit$converged)
            stop("attempt to find suitable starting values failed")
        coefs <- fit$coefficients
        if(any(is.na(coefs))) {
            warning("design appears to be rank-deficient, so dropping some coefs")
            keep <- names(coefs)[!is.na(coefs)]
            coefs <- coefs[keep]
            x <- x[, keep[-1L], drop = FALSE]
            pc <- ncol(x)
        }
        logit <- function(p) log(p/(1 - p))
        spacing <- logit((1L:q)/(q+1L)) # just a guess
        if(method != "logit") spacing <- spacing/1.7
        gammas <- -coefs[1L] + spacing - spacing[q1]
        start <- c(coefs[-1L], gammas)
    } else if(length(start) != pc + q)
	stop("'start' is not of the correct length")

    ans <- newpolr.fit(x, y, wt, start, offset, method, hessian = Hess, ...)
    beta <- ans$coefficients
    zeta <- ans$zeta
    deviance <- ans$deviance
    grad <- ans$grad
    res <- ans$res
    niter <- c(f.evals = res$counts[1L], g.evals = res$counts[2L])

    eta <- if(pc) offset + drop(x %*% beta) else offset + rep(0, n)
    pfun <- switch(method, logit = plogis, probit = pnorm,
                   cloglog = pgumbel, cauchit = pcauchy)
    dfun <- switch(method, logit = dlogis, probit = dnorm,
                   cloglog = dgumbel, cauchit = dcauchy)
    cumpr <- matrix(pfun(matrix(zeta, n, q, byrow=TRUE) - eta), , q)
    dcumpr <- matrix(dfun(matrix(zeta, n, q, byrow=TRUE) - eta), , q)
    dimnames(cumpr) <- dimnames(dcumpr) <- list(row.names(m), names(zeta))
    
    fitted <- t(apply(cumpr, 1L, function(x) diff(c(0, x, 1))))
    dfitted <- t(apply(dcumpr, 1L, function(x) diff(c(0, x, 0))))
    dimnames(fitted) <- dimnames(dfitted) <- list(row.names(m), lev)

    fit <- list(coefficients = beta, zeta = zeta, deviance = deviance,
                grad = grad, cumpr=cumpr, dcumpr=dcumpr,
                fitted.values = fitted, dfitted.values=dfitted,
                lev = lev, terms = Terms,
                df.residual = sum(wt) - pc - q, edf = pc + q, n = sum(wt),
                nobs = sum(wt), call = match.call(), method = method,
		convergence = res$convergence, niter = niter, lp = eta)
    if(Hess) {
        dn <- c(names(beta), names(zeta))
        H <- res$hessian
        dimnames(H) <- list(dn, dn)

        # Shuffle
        Hessian <- matrix(nrow=nrow(H),ncol=ncol(H))
        dimnames(Hessian) <- list(c(names(zeta),names(beta)),c(names(zeta),names(beta)))
        Hessian[names(zeta),names(zeta)] <- -H[names(zeta),names(zeta)]
        Hessian[names(beta),names(beta)] <- -H[names(beta),names(beta)]
        Hessian[names(beta),names(zeta)] <-  H[names(beta),names(zeta)]
        Hessian[names(zeta),names(beta)] <-  H[names(zeta),names(beta)]

        fit$Hessian <- Hessian
    }
    if(model) fit$model <- m
    fit$na.action <- attr(m, "na.action")
    fit$contrasts <- cons
    fit$xlevels <- .getXlevels(Terms, m)
    class(fit) <- "polr"
    fit
}


#' @importFrom stats dlogis dnorm dcauchy optim
newpolr.fit <- function(x, y, wt, start, offset, method, ...)
{
    fmin <- function(beta) {
        theta <- beta[pc + ind_q]
        gamm <- c(-Inf , theta, Inf)
        eta <- offset
        if (pc) eta <- eta + drop(x %*% beta[ind_pc])
        pr <- pfun(pmin(100, gamm[y + 1L] - eta)) -
            pfun(pmax(-100, gamm[y] - eta))
        if (all(pr > 0)) -sum(wt * log(pr)) else Inf
    }

    gmin <- function(beta)
    {
        theta <- beta[pc + ind_q]
        gamm <- c(-Inf, theta, Inf)
        eta <- offset
        if(pc) eta <- eta + drop(x %*% beta[ind_pc])
        
        z1 <- pmin(100, gamm[y+1L] - eta)
        z2 <- pmax(-100, gamm[y] - eta)
        pr <- pfun(z1) - pfun(z2)
        p1 <- dfun(z1); p2 <- dfun(z2)
        g1 <- if(pc) crossprod(x, wt*(p1 - p2)/pr) else numeric()
        xx <- .polrY1*p1 - .polrY2*p2
        g2 <- - crossprod(xx, wt/pr)
        if (all(pr > 0)) c(g1, g2) else rep(NA_real_, pc+q)
    }

    gfun <- function(beta) {
      theta <- beta[pc + ind_q]
      gamm <- c(-Inf, theta, Inf)
      eta <- offset
      if(pc) eta <- eta + drop(x %*% beta[ind_pc])

      z1 <- gamm[y + 1L] - eta
      z2 <- gamm[y] - eta
      pr <- pfun(z1) - pfun(z2)
      p1 <- dfun(z1)
      p2 <- dfun(z2)

      g1 <- if(pc) x * (wt*(p1 - p2)/pr) else numeric()
      xx <- .polrY1*p1 - .polrY2*p2
      g2 <- - xx * (wt/pr)
      if(all(pr > 0)) cbind(-g2, g1) else matrix(NA_real_, n, pc+q)
    }

    dgamma <- function(beta) {
      theta <- beta[pc + ind_q]
      gamm <- c(-Inf, theta, Inf)
      eta <- offset
      if(pc) eta <- eta + drop(x %*% beta[ind_pc])

      p1 <- dfun(gamm[y + 1L] - eta)
      
      g1 <- if(pc) x * (wt*(p1))
      g2 <- .polrY1*wt*p1

      cbind(g1,g2)
    }

    dp0 <- function(beta) {
      theta <- beta[pc + ind_q]
      gamm <- c(-Inf, theta, Inf)
      eta <- offset
      if(pc) eta <- eta + drop(x %*% beta[ind_pc])

      p1 <- dfun(gamm[y + 1L] - eta)
      p2 <- dfun(gamm[y] - eta)
      
      g1 <- if(pc) x * (wt*(p1 - p2))
      g2 <- wt*(.polrY1*p1 - .polrY2*p2)

      cbind(g1,g2)
    }
    
    pfun <- switch(method, logit = plogis, probit = pnorm,
                   cloglog = pgumbel, cauchit = pcauchy)
    dfun <- switch(method, logit = dlogis, probit = dnorm,
                   cloglog = dgumbel, cauchit = dcauchy)
    n <- nrow(x)
    pc <- ncol(x)
    ind_pc <- seq_len(pc)
    lev <- levels(y)
    if(length(lev) <= 2L) stop("response must have 3 or more levels")
    y <- unclass(y)
    q <- length(lev) - 1L
    ind_q <- seq_len(q)
    Y <- matrix(0L, n, q)
    .polrY1 <- col(Y) == y; .polrY2 <- col(Y) == (y - 1L)

    ## pc could be 0.
    res <- optim(start, fmin, gmin, method="BFGS", ...)
    beta <- res$par[ind_pc]
    theta <- res$par[pc + ind_q]
    zeta <- theta

    grad <- gfun(res$par)
    
    deviance <- 2 * res$value
    names(zeta) <- paste(lev[-length(lev)], lev[-1L], sep="|")
    if(pc) names(beta) <- colnames(x)

    dn <- c(names(beta), names(zeta))
    colnames(grad) <- dn
    
    list(coefficients = beta, zeta = zeta, grad=grad, deviance = deviance,
         res = res)
}
