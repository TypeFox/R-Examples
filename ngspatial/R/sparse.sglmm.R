
sparse.sglmm.fit.binomial = function(Y, X, A, M, family, beta.start, V, offset, tol, minit, maxit, sigma.s, sigma.b, verbose)
{
    iterations = minit
    n = length(Y)
    Q = t(M) %*% (diag(rowSums(A), n) - A) %*% M
    p = ncol(X)
    beta = matrix(0, iterations + 1, p)
    beta[1, ] = beta.start
    tau.s = numeric(iterations + 1)
    tau.s[1] = 0.1
    q = ncol(M)
    gamma = matrix(0, iterations + 1, q)
    gamma[1, ] = rnorm(q, 0, 1)
    a.s = 0.5
    b.s = 2000
    linkinv = family$linkinv
    if (is.null(offset))
        offset = rep(0, n)
    start = 1
    k = 1
    five.pct = round(0.05 * maxit, 0)
    repeat
    {
        for (j in (start + 1):(start + iterations))
        {
            k = k + 1
            if (verbose && k %% five.pct == 0)
                cat("Progress => ", round(k / maxit * 100, 0), "%\n", sep = "")
            b = beta[j - 1, ]
            b_ = V %*% rnorm(p) + b
            gam = gamma[j - 1, ]
            eta_ = offset + X %*% b_ + M %*% gam
            eta = offset + X %*% b + M %*% gam
            logAlpha = t(Y) %*% (log(linkinv(eta_)) - log(linkinv(eta)))
            logAlpha = logAlpha + t(1 - Y) %*% (log(1 - linkinv(eta_)) - log(1 - linkinv(eta)))
            logAlpha = logAlpha + (sum(b^2) - sum(b_^2)) / (2 * sigma.b^2)
            if (log(runif(1)) < logAlpha)
                b = b_
            beta[j, ] = b
            gam_ = rnorm(q, 0, sigma.s) + gam
            eta_ = offset + X %*% b + M %*% gam_
            eta = offset + X %*% b + M %*% gam
            logAlpha = t(Y) %*% (log(linkinv(eta_)) - log(linkinv(eta)))
            logAlpha = logAlpha + t(1 - Y) %*% (log(1 - linkinv(eta_)) - log(1 - linkinv(eta)))
            logAlpha = logAlpha + 0.5 * tau.s[j - 1] * (t(gam) %*% Q %*% gam - t(gam_) %*% Q %*% gam_)
            if (log(runif(1)) < logAlpha)
                gam = gam_
            gamma[j, ] = gam
            tau.s[j] = rgamma(1, a.s + q / 2, t(gam) %*% Q %*% gam / 2 + 1 / b.s)
            if (j == maxit)
                break
        }
        if (j == maxit)
        {
            beta = as.matrix(beta[1:maxit, ])
            gamma = as.matrix(gamma[1:maxit, ])
            tau.s = tau.s[1:maxit]
            break
        }
        done = TRUE
        for (j in 1:p)
        {
            temp = bm(beta[, j])
            if (temp$se > tol)
            {
                done = FALSE
                break
            }
        }
        if (done)
        {
            for (j in 1:q)
            {
                temp = bm(gamma[, j])
                if (temp$se > tol)
                {
                    done = FALSE
                    break
                }
            }
        }
        if (done)
        {
            temp = bm(tau.s)
            if (temp$se > tol)
                done = FALSE
        }
        if (done)
            break
        else
        {
            start = start + iterations
            temp = matrix(0, iterations, p)
            beta = rbind(beta, temp)
            temp = matrix(0, iterations, q)
            gamma = rbind(gamma, temp)
            tau.s = c(tau.s, rep(0, iterations))
        }
    }
    coefficients = numeric(p)
    beta.mcse = numeric(p)
    names(coefficients) = names(beta.mcse) = colnames(X)
    for (j in 1:p)
    {
        temp = bm(beta[, j])
        coefficients[j] = temp$est
        beta.mcse[j] = temp$se
    }
    gamma.est = numeric(q)
    gamma.mcse = numeric(q)
    for (j in 1:q)
    {
        temp = bm(gamma[, j])
        gamma.est[j] = temp$est
        gamma.mcse[j] = temp$se
    }
    temp = bm(tau.s)
    tau.s.est = temp$est
    tau.s.mcse = temp$se
    linear.predictors = numeric(n)
    fitted.values = numeric(n)
    iter = length(tau.s)
    v = numeric(iter)
    for (j in 1:iter)
    {
        eta = offset + X %*% beta[j, ] + M %*% gamma[j, ]
        linear.predictors = linear.predictors + eta / iter
        p = linkinv(eta)
        fitted.values = fitted.values + p / iter
        v[j] = -2 * sum(dbinom(Y, 1, p, log = TRUE))
    }
    D.bar = mean(v)
    pD = D.bar + 2 * sum(dbinom(Y, 1, fitted.values, log = TRUE))
    dic = D.bar + pD
    residuals = Y - fitted.values
    beta.accept = sum(diff(beta[, 1]) != 0) / iter
    gamma.accept = sum(diff(gamma[, 1]) != 0) / iter
    object = list(coefficients = coefficients, fitted.values = fitted.values,
                  linear.predictors = linear.predictors, residuals = residuals,
                  beta.sample = beta, gamma.sample = gamma, tau.s.sample = tau.s,
                  beta.mcse = beta.mcse, gamma.mcse = gamma.mcse, tau.s.mcse = tau.s.mcse,
                  gamma.est = gamma.est, tau.s.est = tau.s.est, iter = iter, dic = dic,
                  D.bar = D.bar, pD = pD, beta.accept = beta.accept, gamma.accept = gamma.accept)
    class(object) = c("sparse.sglmm")
    object
}

sparse.sglmm.fit.gaussian = function(Y, X, A, M, beta.start, offset, tol, minit, maxit, hyper, verbose)
{
    iterations = minit
    n = length(Y)
    Q = t(M) %*% (diag(rowSums(A), n) - A) %*% M
    p = ncol(X)
    XtX = t(X) %*% X
    MtM = t(M) %*% M
    beta = matrix(0, iterations + 1, p)
    beta[1, ] = beta.start
    tau.s = numeric(iterations + 1)
    tau.s[1] = 0.1
    tau.h = numeric(iterations + 1)
    tau.h[1] = 1
    q = ncol(M)
    gamma = matrix(0, iterations + 1, q)
    gamma[1, ] = rnorm(q, 0, 1)
    a.s = 0.5
    b.s = 2000
    a.h = hyper$a.h
    b.h = hyper$b.h
    if (is.null(offset))
        offset = rep(0, n)
    start = 1
    k = 1
    five.pct = round(0.05 * maxit, 0)
    K = diag(1 / hyper$sigma.b^2, p)
    repeat
    {
        for (j in (start + 1):(start + iterations))
        {
            k = k + 1
            if (verbose && k %% five.pct == 0)
                cat("Progress => ", round(k / maxit * 100, 0), "%\n", sep = "")
            V = solve(K + tau.h[j - 1] * XtX)
            mu = V %*% (tau.h[j - 1] * t(X) %*% (Y - M %*% gamma[j - 1, ]))
            beta[j, ] = t(chol(V)) %*% rnorm(p) + mu
            V = solve(tau.s[j - 1] * Q + tau.h[j - 1] * MtM)
            mu = V %*% (tau.h[j - 1] * t(M) %*% (Y - X %*% beta[j, ]))
            gamma[j, ] = t(chol(V)) %*% rnorm(q) + mu
            b = 0.5 * sum((Y - X %*% beta[j, ] - M %*% gamma[j, ])^2) + 1 / b.h
            tau.h[j] = rgamma(1, a.h + n / 2, b)
            b = 0.5 * t(gamma[j, ]) %*% Q %*% gamma[j, ] + 1 / b.s
            tau.s[j] = rgamma(1, a.s + q / 2, b)
            if (j == maxit)
                break
        }
        if (j == maxit)
        {
            beta = as.matrix(beta[1:maxit, ])
            gamma = as.matrix(gamma[1:maxit, ])
            tau.s = tau.s[1:maxit]
            tau.h = tau.h[1:maxit]
            break
        }
        done = TRUE
        for (j in 1:p)
        {
            temp = bm(beta[, j])
            if (temp$se > tol)
            {
                done = FALSE
                break
            }
        }
        if (done)
        {
            for (j in 1:q)
            {
                temp = bm(gamma[, j])
                if (temp$se > tol)
                {
                    done = FALSE
                    break
                }
            }
        }
        if (done)
        {
            temp = bm(tau.s)
            if (temp$se > tol)
                done = FALSE
        }
        if (done)
        {
            temp = bm(tau.h)
            if (temp$se > tol)
                done = FALSE
        }
        if (done)
            break
        else
        {
            start = start + iterations
            temp = matrix(0, iterations, p)
            beta = rbind(beta, temp)
            temp = matrix(0, iterations, q)
            gamma = rbind(gamma, temp)
            tau.s = c(tau.s, rep(0, iterations))
            tau.h = c(tau.h, rep(0, iterations))
        }
    }
    coefficients = numeric(p)
    beta.mcse = numeric(p)
    names(coefficients) = names(beta.mcse) = colnames(X)
    for (j in 1:p)
    {
        temp = bm(beta[, j])
        coefficients[j] = temp$est
        beta.mcse[j] = temp$se
    }
    gamma.est = numeric(q)
    gamma.mcse = numeric(q)
    for (j in 1:q)
    {
        temp = bm(gamma[, j])
        gamma.est[j] = temp$est
        gamma.mcse[j] = temp$se
    }
    temp = bm(tau.s)
    tau.s.est = temp$est
    tau.s.mcse = temp$se
    temp = bm(tau.h)
    tau.h.est = temp$est
    tau.h.mcse = temp$se
    linear.predictors = numeric(n)
    iter = length(tau.s)
    v = numeric(iter)
    for (j in 1:iter)
    {
        mu = offset + X %*% beta[j, ] + M %*% gamma[j, ]
        linear.predictors = linear.predictors + mu / iter
        v[j] = -2 * sum(dnorm(Y, mu, 1 / sqrt(tau.h[j]), log = TRUE))
    }
    D.bar = mean(v)
    pD = D.bar + 2 * sum(dnorm(Y, linear.predictors, 1 / sqrt(tau.h.est), log = TRUE))
    dic = D.bar + pD
    fitted.values = linear.predictors
    residuals = Y - fitted.values
    object = list(coefficients = coefficients, fitted.values = fitted.values,
                  linear.predictors = linear.predictors, residuals = residuals,
                  beta.sample = beta, gamma.sample = gamma, tau.s.sample = tau.s,
                  tau.h.sample = tau.h, tau.h.mcse = tau.h.mcse, beta.mcse = beta.mcse,
                  gamma.mcse = gamma.mcse, tau.s.mcse = tau.s.mcse, gamma.est = gamma.est,
                  tau.s.est = tau.s.est, iter = iter, dic = dic, tau.h.est = tau.h.est,
                  D.bar = D.bar, pD = pD)
    class(object) = c("sparse.sglmm")
    object
}

sparse.sglmm.fit.poisson = function(Y, X, A, M, family, beta.start, V, offset, tol, minit, maxit, tune, hyper, verbose)
{
    iterations = minit
    n = length(Y)
    Q = t(M) %*% (diag(rowSums(A), n) - A) %*% M
    p = ncol(X)
    beta = matrix(0, iterations + 1, p)
    beta[1, ] = beta.start
    tau.s = numeric(iterations + 1)
    tau.s[1] = 0.1
    q = ncol(M)
    gamma = matrix(0, iterations + 1, q)
    gamma[1, ] = rnorm(q, 0, 1)
    sigma.s = tune$sigma.s
    a.s = 0.5
    b.s = 2000
    hetero = FALSE
    if (! is.null(tune$sigma.h))
    {
        R = t(M) %*% M
        delta = matrix(0, iterations + 1, q)
        delta[1, ] = rnorm(q, 0, 1)
        sigma.h = tune$sigma.h
        a.h = hyper$a.h
        b.h = hyper$b.h
        hetero = TRUE
        tau.h = numeric(iterations + 1)
        tau.h[1] = 1
    }
    sigma.b = hyper$sigma.b
    linkinv = family$linkinv
    if (is.null(offset))
        offset = rep(0, n)
    start = 1
    k = 1
    five.pct = round(0.05 * maxit, 0)
    repeat
    {
        if (hetero)
        {
            for (j in (start + 1):(start + iterations))
            {
                k = k + 1
                if (verbose && k %% five.pct == 0)
                    cat("Progress => ", round(k / maxit * 100, 0), "%\n", sep = "")
                b = beta[j - 1, ]
                b_ = V %*% rnorm(p) + b
                gam = gamma[j - 1, ]
                delt = delta[j - 1, ]
                eta_ = offset + X %*% b_ + M %*% gam + M %*% delt
                eta = offset + X %*% b + M %*% gam + M %*% delt
                logAlpha = t(Y) %*% (log(linkinv(eta_)) - log(linkinv(eta)))
                logAlpha = logAlpha + sum(linkinv(eta)) - sum(linkinv(eta_))
                logAlpha = logAlpha + (sum(b^2) - sum(b_^2)) / (2 * sigma.b^2)
                if (log(runif(1)) < logAlpha)
                    b = b_
                beta[j, ] = b
                gam_ = rnorm(q, 0, sigma.s) + gam
                eta_ = offset + X %*% b + M %*% gam_ + M %*% delt
                eta = offset + X %*% b + M %*% gam + M %*% delt
                logAlpha = t(Y) %*% (log(linkinv(eta_)) - log(linkinv(eta)))
                logAlpha = logAlpha + sum(linkinv(eta)) - sum(linkinv(eta_))
                logAlpha = logAlpha + 0.5 * tau.s[j - 1] * (t(gam) %*% Q %*% gam - t(gam_) %*% Q %*% gam_)
                if (log(runif(1)) < logAlpha)
                    gam = gam_
                gamma[j, ] = gam
                tau.s[j] = rgamma(1, a.s + q / 2, t(gam) %*% Q %*% gam / 2 + 1 / b.s)
                delt_ = rnorm(q, 0, sigma.h) + delt
                eta_ = offset + X %*% b + M %*% gam + M %*% delt_
                eta = offset + X %*% b + M %*% gam + M %*% delt
                logAlpha = t(Y) %*% (log(linkinv(eta_)) - log(linkinv(eta)))
                logAlpha = logAlpha + sum(linkinv(eta)) - sum(linkinv(eta_))
                logAlpha = logAlpha + 0.5 * tau.h[j - 1] * (t(delt) %*% R %*% delt - t(delt_) %*% R %*% delt_)
                if (log(runif(1)) < logAlpha)
                    delt = delt_
                delta[j, ] = delt
                tau.h[j] = rgamma(1, a.h + q / 2, t(delt) %*% R %*% delt / 2 + 1 / b.h)
                if (j == maxit)
                    break
            }
        }
        else
        {
            for (j in (start + 1):(start + iterations))
            {
                k = k + 1
                if (verbose && k %% five.pct == 0)
                    cat("Progress => ", round(k / maxit * 100, 0), "%\n", sep = "")
                b = beta[j - 1, ]
                b_ = V %*% rnorm(p) + b
                gam = gamma[j - 1, ]
                eta_ = offset + X %*% b_ + M %*% gam
                eta = offset + X %*% b + M %*% gam
                logAlpha = t(Y) %*% (log(linkinv(eta_)) - log(linkinv(eta)))
                logAlpha = logAlpha + sum(linkinv(eta)) - sum(linkinv(eta_))
                logAlpha = logAlpha + (sum(b^2) - sum(b_^2)) / 2000000
                if (log(runif(1)) < logAlpha)
                    b = b_
                beta[j, ] = b
                gam_ = rnorm(q, 0, sigma.s) + gam
                eta_ = offset + X %*% b + M %*% gam_
                eta = offset + X %*% b + M %*% gam
                logAlpha = t(Y) %*% (log(linkinv(eta_)) - log(linkinv(eta)))
                logAlpha = logAlpha + sum(linkinv(eta)) - sum(linkinv(eta_))
                logAlpha = logAlpha + 0.5 * tau.s[j - 1] * (t(gam) %*% Q %*% gam - t(gam_) %*% Q %*% gam_)
                if (log(runif(1)) < logAlpha)
                    gam = gam_
                gamma[j, ] = gam
                tau.s[j] = rgamma(1, a.s + q / 2, t(gam) %*% Q %*% gam / 2 + 1 / b.s)
                if (j == maxit)
                    break
            }

        }
        if (j == maxit)
        {
            beta = as.matrix(beta[1:maxit, ])
            gamma = as.matrix(gamma[1:maxit, ])
            tau.s = tau.s[1:maxit]
            if (hetero)
            {
                delta = as.matrix(delta[1:maxit, ])
                tau.h = tau.h[1:maxit]
            }
            break
        }
        done = TRUE
        for (j in 1:p)
        {
            temp = bm(beta[, j])
            if (temp$se > tol)
            {
                done = FALSE
                break
            }
        }
        if (done)
        {
            for (j in 1:q)
            {
                temp = bm(gamma[, j])
                if (temp$se > tol)
                {
                    done = FALSE
                    break
                }
            }
        }
        if (done)
        {
            temp = bm(tau.s)
            if (temp$se > tol)
                done = FALSE
        }
        if (done && hetero)
        {
            temp = bm(tau.h)
            if (temp$se > tol)
                done = FALSE
        }
        if (done)
            break
        else
        {
            start = start + iterations
            temp = matrix(0, iterations, p)
            beta = rbind(beta, temp)
            temp = matrix(0, iterations, q)
            gamma = rbind(gamma, temp)
            tau.s = c(tau.s, rep(0, iterations))
            if (hetero)
            {
                delta = rbind(delta, temp)
                tau.h = c(tau.h, rep(0, iterations))
            }
        }
    }
    coefficients = numeric(p)
    beta.mcse = numeric(p)
    names(coefficients) = names(beta.mcse) = colnames(X)
    for (j in 1:p)
    {
        temp = bm(beta[, j])
        coefficients[j] = temp$est
        beta.mcse[j] = temp$se
    }
    gamma.est = numeric(q)
    gamma.mcse = numeric(q)
    for (j in 1:q)
    {
        temp = bm(gamma[, j])
        gamma.est[j] = temp$est
        gamma.mcse[j] = temp$se
    }
    temp = bm(tau.s)
    tau.s.est = temp$est
    tau.s.mcse = temp$se
    if (hetero)
    {
        delta.est = numeric(q)
        delta.mcse = numeric(q)
        for (j in 1:q)
        {
            temp = bm(delta[, j])
            delta.est[j] = temp$est
            delta.mcse[j] = temp$se
        }
        temp = bm(tau.h)
        tau.h.est = temp$est
        tau.h.mcse = temp$se
    }
    linear.predictors = numeric(n)
    fitted.values = numeric(n)
    iter = length(tau.s)
    v = numeric(iter)
    for (j in 1:iter)
    {
        eta = offset + X %*% beta[j, ] + M %*% gamma[j, ]
        if (hetero)
            eta = eta + M %*% delta[j, ]
        linear.predictors = linear.predictors + eta / iter
        lambda = linkinv(eta)
        fitted.values = fitted.values + lambda / iter
        v[j] = -2 * sum(dpois(Y, lambda, log = TRUE))
    }
    D.bar = mean(v)
    pD = D.bar + 2 * sum(dpois(Y, fitted.values, log = TRUE))
    dic = D.bar + pD
    residuals = Y - fitted.values
    beta.accept = sum(diff(beta[, 1]) != 0) / iter
    gamma.accept = sum(diff(gamma[, 1]) != 0) / iter
    if (hetero)
        delta.accept = sum(diff(delta[, 1]) != 0) / iter
    object = list(coefficients = coefficients, fitted.values = fitted.values,
                  linear.predictors = linear.predictors, residuals = residuals,
                  beta.sample = beta, gamma.sample = gamma, tau.s.sample = tau.s,
                  beta.mcse = beta.mcse, gamma.mcse = gamma.mcse, tau.s.mcse = tau.s.mcse,
                  gamma.est = gamma.est, tau.s.est = tau.s.est, iter = iter, D.bar = D.bar,
                  pD = pD, dic = dic, beta.accept = beta.accept, gamma.accept = gamma.accept)
    if (hetero)
        object = c(object,
                 list(delta.sample = delta, tau.h.sample = tau.h, delta.mcse = delta.mcse,
                      tau.h.mcse = tau.h.mcse, delta.est = delta.est, tau.h.est = tau.h.est,
                      delta.accept = delta.accept))
    class(object) = c("sparse.sglmm")
    object
}

#' Fit a sparse SGLMM.
#'
#' @details This function fits the sparse areal SGLMM of Hughes and Haran (2013). The first stage of the model is \deqn{g(\mu_i)=x_i^\prime\beta+m_i^\prime\gamma\hspace{1 cm}(i=1,\dots,n)}{g(\mu_i)=x_i'\beta+m_i'\gamma   (i=1,\dots,n)} or, in vectorized form, \deqn{g(\mu)=X\beta+M\gamma,} where \eqn{X} is the design matrix, \eqn{\beta} is a \eqn{p}-vector of regression coefficients, the columns of \eqn{M} are \eqn{q} eigenvectors of the Moran operator, and \eqn{\gamma} are spatial random effects. Arguments \code{attractive} and \code{repulsive} can be used to control the number of eigenvectors used. The default values are 50 and 0, respectively, which corresponds to pure spatial smoothing. Inclusion of some repulsive eigenvectors can be advantageous in certain applications.\cr\cr The second stage, i.e., the prior for \eqn{\gamma}, is \deqn{p(\gamma\mid\tau_s)\propto\tau_s^{q/2}\exp\left(-\frac{\tau_s}{2}\gamma^\prime M^\prime QM\gamma\right),}{p(\gamma | \tau_s) proportional to \tau_s^(q/2)exp(-\tau_s/2 \gamma'M'QM\gamma'),} where \eqn{\tau_s} is a smoothing parameter and \eqn{Q} is the graph Laplacian.\cr\cr The prior for \eqn{\beta} is spherical \eqn{p}-variate normal with mean zero and common standard deviation \code{sigma.b}, which defaults to 1,000. The prior for \eqn{\tau_s} is gamma with parameters 0.5 and 2,000.\cr\cr When the response is normally distributed, the identity link is assumed, in which case the first stage is \deqn{\mu=X\beta+M\gamma+M\delta,} where \eqn{\delta} are heterogeneity random effects. When the response is Poisson distributed, heterogeneity random effects are optional. In any case, the prior on \eqn{\delta} is spherical \eqn{q}-variate normal with mean zero and common variance \eqn{1/\tau_h}. The prior for \eqn{\tau_h} is gamma with parameters \eqn{a_h} and \eqn{b_h}, the values of which are controlled by the user through argument \code{hyper}.\cr\cr If the response is Bernoulli or Poisson, \eqn{\beta} and \eqn{\gamma} are updated using Metropolis-Hastings random walks with normal proposals. The proposal covariance matrix for \eqn{\beta} is the estimated asymptotic covariance matrix from a \code{\link{glm}} fit to the data (see \code{\link{vcov}}). The proposal for \eqn{\gamma} is spherical normal with common standard deviation \code{sigma.s}.\cr\cr The updates for \eqn{\tau_s} and \eqn{\tau_h} are Gibbs updates irrespective of the response distribution.\cr\cr If the response is Poisson distributed and heterogeneity random effects are included, those random effects are updated using a Metropolis-Hastings random walk with a spherical normal proposal. The common standard deviation is \code{sigma.h}.\cr\cr If the response is normally distributed, all updates are Gibbs updates.
#' @param formula an object of class \code{\link{formula}}: a symbolic description of the model to be fitted.
#' @param family a description of the error distribution and link function to be used in the model. This can be a character string naming a family function, a family function, or the result of a call to a family function. (See \code{\link{family}} for details of family functions.) Supported families are \code{gaussian} (default), \code{binomial}, and \code{poisson}.
#' @param data an optional data frame, list, or environment (or object coercible by \code{\link{as.data.frame}} to a data frame) containing the variables in the model. If not found in \code{data}, the variables are taken from \code{environment(formula)}, typically the environment from which \code{sparse.sglmm} is called.
#' @param offset this can be used to specify an \emph{a priori} known component to be included in the linear predictor during fitting. This should be \code{NULL} or a numeric vector of length equal to the number of cases. One or more \code{\link{offset}} terms can be included in the formula instead or as well, and if more than one is specified their sum is used. See \code{\link{model.offset}}.
#' @param A the adjacency matrix for the underlying graph.
#' @param attractive the number of attractive Moran eigenvectors to use. The default is 50. See `Details' for more information.
#' @param repulsive the number of repulsive Moran eigenvectors to use. The default is 0. See `Details' for more information.
#' @param tol a tolerance. If all Monte Carlo standard errors are smaller than \code{tol}, no more samples are drawn from the posterior. The default is 0.01.
#' @param minit the minimum sample size. This should be large enough to permit accurate estimation of Monte Carlo standard errors. The default is 10,000.
#' @param maxit the maximum sample size. Sampling from the posterior terminates when all Monte Carlo standard errors are smaller than \code{tol} or when \code{maxit} samples have been drawn, whichever happens first. The default is 1,000,000.
#' @param tune (where relevant) a list containing \code{sigma.s} and \code{sigma.h}. These are the standard deviations for the \eqn{\gamma} and \eqn{\delta} proposals, respectively.
#' @param hyper a list containing \code{sigma.b}, the prior standard deviation for \eqn{\beta}, and (where relevant)  \code{a.h} and \code{b.h}, the parameters of the gamma prior for \eqn{\tau_h}.
#' @param model a logical value indicating whether the model frame should be included as a component of the returned value.
#' @param x a logical value indicating whether the model matrix used in the fitting process should be returned as a component of the returned value.
#' @param y a logical value indicating whether the response vector used in the fitting process should be returned as a component of the returned value.
#' @param verbose a logical value indicating whether to print MCMC progress to the screen. Defaults to \code{FALSE}.
#' @return \code{sparse.sglmm} returns an object of class \dQuote{\code{sparse.sglmm}}, which is a list containing the following components.
#'         \item{coefficients}{the estimated regression coefficients.}
#'         \item{fitted.values}{the fitted mean values, obtained by transforming the linear predictors by the inverse of the link function.}
#'         \item{linear.predictors}{the linear fit on link scale.}
#'         \item{residuals}{the response residuals.}
#'         \item{iter}{the size of the posterior sample.}
#'         \item{beta.sample}{an \code{iter} by \eqn{p} matrix containing the posterior samples for \eqn{\beta}.}
#'         \item{gamma.sample}{an \code{iter} by \eqn{q} matrix containing the posterior samples for \eqn{\gamma}.}
#'         \item{delta.sample}{(where relevant) an \code{iter} by \eqn{q} matrix containing the posterior samples for \eqn{\delta}.}
#'         \item{tau.s.sample}{a vector containing the posterior samples for \eqn{\tau_s}.}
#'         \item{tau.h.sample}{(where relevant) a vector containing the posterior samples for \eqn{\tau_h}.}
#'         \item{gamma.est}{the estimate of \eqn{\gamma}.}
#'         \item{delta.est}{(where relevant) the estimate of \eqn{\delta}.}
#'         \item{tau.s.est}{the estimate of \eqn{\tau_s}.}
#'         \item{tau.h.est}{(where relevant) the estimate of \eqn{\tau_h}.}
#'         \item{beta.mcse}{the Monte Carlo standard errors for \eqn{\beta}.}
#'         \item{gamma.mcse}{the Monte Carlo standard errors for \eqn{\gamma}.}
#'         \item{delta.mcse}{(where relevant) the Monte Carlo standard errors for \eqn{\delta}.}
#'         \item{tau.s.mcse}{the Monte Carlo standard error for \eqn{\tau_s}.}
#'         \item{tau.h.mcse}{(where relevant) the Monte Carlo standard error for \eqn{\tau_h}.}
#'         \item{D.bar}{the goodness of fit component of the DIC.}
#'         \item{pD}{the penalty component of the DIC.}
#'         \item{dic}{the deviance information criterion.}
#'         \item{beta.accept}{the acceptance rate for \eqn{\beta}.}
#'         \item{gamma.accept}{the acceptance rate for \eqn{\gamma}.}
#'         \item{delta.accept}{(where relevant) the acceptance rate for \eqn{\delta}.}
#'         \item{y}{if requested (the default), the \code{y} vector used.}
#'         \item{X}{if requested, the model matrix.}
#'         \item{M}{if requested, the matrix of Moran eigenvectors.}
#'         \item{eigen.values}{if requested, the spectrum of the Moran operator.}
#'         \item{hyper}{a list containing the names and values of the hyperparameters.}
#'         \item{tune}{a list containing the names and values of the tuning parameters.}
#'         \item{model}{if requested (the default), the model frame.}
#'         \item{call}{the matched call.}
#'         \item{formula}{the formula supplied.}
#'         \item{terms}{the \code{\link{terms}} object used.}
#'         \item{data}{the \code{data} argument.}
#'         \item{offset}{the offset vector used.}
#'         \item{xlevels}{(where relevant) a record of the levels of the factors used in fitting.}
#' @references
#' Hughes, J. and Haran, M. (2013) Dimension reduction and alleviation of confounding for spatial generalized linear mixed models. \emph{Journal of the Royal Statistical Society, Series B}, \bold{75}(1), 139--159.
#' @seealso \code{\link{residuals.sparse.sglmm}}, \code{\link{summary.sparse.sglmm}}, \code{\link{vcov.sparse.sglmm}}
#' @export
#' @examples \dontrun{
#'
#' The following code duplicates the analysis described in (Hughes and Haran, 2013). The data are
#' infant mortality data for 3,071 US counties. We do a spatial Poisson regression with offset.
#'
#' data(infant)
#' infant$low_weight = infant$low_weight / infant$births
#' attach(infant)
#' Z = deaths
#' X = cbind(1, low_weight, black, hispanic, gini, affluence, stability)
#' data(A)
#' set.seed(123456)
#' fit = sparse.sglmm(Z ~ X - 1 + offset(log(births)), family = poisson, A = A,
#'                    tune = list(sigma.s = 0.02), verbose = TRUE)
#' summary(fit)
#' } 

sparse.sglmm = function(formula, family = gaussian, data, offset, A, attractive = 50, repulsive = 0, tol = 0.01, minit = 10000,
	                    maxit = 1000000, tune = list(), hyper = list(), model = TRUE, x = FALSE, y = FALSE, verbose = FALSE)
{
    cl = match.call()
    if (is.character(family)) 
        family = get(family, mode = "function", envir = parent.frame())
    if (is.function(family)) 
        family = family()
    if (is.null(family$family))
    {
        print(family)
        stop("'family' was not recognized.")
    }
    if (missing(data))
        data = environment(formula)
    nonspat = if (missing(offset)) glm(formula, family, data)
              else glm(formula, family, data, offset = offset)
    wts = 1 / sqrt(weights(nonspat, type = "working"))
    if (! family$family %in% c("binomial", "gaussian", "poisson"))
        stop("'family' must be binomial, gaussian, or poisson.")
    if (family$family == "gaussian" && family$link != "identity")
        stop("For the gaussian family, only the identity link is supported.")
    if (missing(A) || ! is.matrix(A) || ! isSymmetric(A) || ! (A == 0 || A == 1))
        stop("You must supply a symmetric binary adjacency matrix.")
    diag(A) = 0
    if (! is.numeric(tol) || length(tol) > 1 || tol <= 0 || tol >= 1)
        stop("'tol' must be a number between 0 and 1.")
    if (! is.numeric(minit) || length(minit) > 1 || ! is.wholenumber(minit) || minit < 1)
        stop("'minit' must be a positive whole number.")
    if (! is.numeric(maxit) || length(maxit) > 1 || ! is.wholenumber(maxit) || maxit < minit)
        stop("\n'maxit' must be a positive whole number greater than or equal to 'minit'.")
    if (! is.list(tune))
        stop("'tune' must be a list.")
    if (family$family != "gaussian")
    {
    	sigma.s = tune$sigma.s
    	if (is.null(sigma.s) || ! is.numeric(sigma.s) || length(sigma.s) > 1 || sigma.s <= 0)
        {
		    if (verbose)
                cat("\nTuning parameter 'sigma.s' must be a positive number. Setting it to the default value of 0.01.\n")
            sigma.s = 0.01
        }
    }
    tune.new = list()
    if (family$family != "gaussian")
        tune.new$sigma.s = sigma.s
    if (! is.list(hyper))
        stop("'hyper' must be a list.")
    sigma.b = hyper$sigma.b
    if (is.null(sigma.b) || ! is.numeric(sigma.b) || length(sigma.b) > 1 || sigma.b <= 0)
    {
	    if (verbose)
            cat("\nHyperparameter 'sigma.b' must be a positive number. Setting it to the default value of 1,000.\n")
        sigma.b = 1000
    }
    hyper.new = list()
    hyper.new$sigma.b = sigma.b
    sigma.h = tune$sigma.h
    if (! is.null(sigma.h) && (! is.numeric(sigma.h) || length(sigma.h) > 1 || sigma.h <= 0))
    {
        if (verbose)
            cat("\nTuning parameter 'sigma.h' must be a positive number. Setting it to the default value of 0.01.\n")
        sigma.h = 0.01
    }
    a.h = hyper$a.h
    b.h = hyper$b.h
    if (family$family == "gaussian" || (family$family == "poisson" && ! is.null(tune$sigma.h)))
    {
	     if (is.null(a.h) || ! is.numeric(a.h) || length(a.h) > 1 || a.h <= 0)
	     {
		     if (verbose)
	             cat("\nHyperparameter 'a.h' must be a positive number. Setting it to the default value of 0.01.\n")
	         a.h = 0.01
	     }
	     if (is.null(b.h) || ! is.numeric(b.h) || length(b.h) > 1 || b.h <= 0)
	     {
		     if (verbose)
	             cat("\nHyperparameter 'b.h' must be a positive number. Setting it to the default value of 100.\n")
	         b.h = 100
	     }
    }
    if (family$family == "poisson" && ! is.null(sigma.h))
        tune.new$sigma.h = sigma.h
    if (family$family == "gaussian" || (family$family == "poisson" && ! is.null(sigma.h)))
    {
    	hyper.new$a.h = a.h
    	hyper.new$b.h = b.h
    }
    tune = tune.new
    hyper = hyper.new
    mf = match.call(expand.dots = FALSE)
    m = match(c("formula", "data", "offset"), names(mf), 0)
    mf = mf[c(1, m)]
    mf[[1]] = as.name("model.frame")
    mf = eval(mf, parent.frame())
    mt = attr(mf, "terms")
    Y = model.response(mf, "numeric")
    n = length(Y)
    offset = as.vector(model.offset(mf))
    X = diag(wts, n) %*% model.matrix(mt, mf)
    P = diag(1, n) - X %*% solve(t(X) %*% X) %*% t(X)
    X = model.matrix(mt, mf)
    if (sum(c(nrow(X), nrow(A)) != length(Y)) > 0)
        stop("The supplied response vector/design matrix/adjacency matrix are not conformable.")
    if (! is.logical(verbose) || length(verbose) > 1)
        stop("'verbose' must be a logical value.")
    if (verbose)
    {
        cat("\nThe Moran operator is being eigendecomposed. This operation may be time consuming.\n")
        flush.console()
        Sys.sleep(0.5)
    }
    eig = eigen(P %*% A %*% P, symmetric = TRUE)
    eigenvalues = eig$values
    maxatt = match(TRUE, sapply(eigenvalues, is.zero)) - 1
    if (! is.numeric(attractive) || length(attractive) > 1 || ! is.wholenumber(attractive) || attractive < 1 || attractive > maxatt)
        stop(gettextf("'attractive' is %d but should be a whole number between 1 and %d.", attractive, maxatt), domain = NA)
    maxrep = match(TRUE, sapply(rev(eigenvalues), is.zero)) - 1
    if (! is.numeric(repulsive) || length(repulsive) > 1 || ! is.wholenumber(repulsive) || repulsive < 0 || repulsive > maxrep)
        stop(gettextf("'repulsive' is %d but should be a whole number between 0 and %d.", repulsive, maxrep), domain = NA)
    M = eig$vectors[, 1:attractive]
	if (repulsive > 0)
	{
		pos2 = length(eigenvalues)
		pos1 = pos2 - repulsive + 1
		M = cbind(M, eig$vectors[, pos1:pos2])
	}
    rm(eig)
    V = t(chol(vcov(nonspat)))
    if (verbose)
    {
        cat("\nWarning: MCMC may be time consuming.\n\n")
        flush.console()
        Sys.sleep(0.5)
    }
    fit = if (family$family == "binomial")
              sparse.sglmm.fit.binomial(Y, X, A, M, family, nonspat$coef, V, offset, tol, minit, maxit, tune$sigma.s, hyper$sigma.b, verbose)
          else if (family$family == "gaussian")
              sparse.sglmm.fit.gaussian(Y, X, A, M, nonspat$coef, offset, tol, minit, maxit, hyper, verbose)
          else sparse.sglmm.fit.poisson(Y, X, A, M, family, nonspat$coef, V, offset, tol, minit, maxit, tune, hyper, verbose)
    fit$xlevels = .getXlevels(mt, mf)
    fit$call = cl
    fit$terms = mt
    fit$formula = formula
    fit$family = family
    fit$offset = offset
    if (model)
        fit$model = mf
    if (x)
    {
        fit$X = X
        fit$M = M
        fit$eigenvalues = eigenvalues
    }
    if (y)
        fit$y = Y
    fit$tune = tune
    fit$hyper = hyper
    fit
}

#' Return the covariance matrix of the regression parameters of a \code{sparse.sglmm} model object.
#'
#' @param object a fitted \code{sparse.sglmm} model object.
#' @param \dots additional arguments.
#' @return An estimate of the posterior covariance matrix of the regression coefficients.
#' @method vcov sparse.sglmm
#' @export

vcov.sparse.sglmm = function(object, ...)
{
    V = cov(object$beta.sample)
    rownames(V) = colnames(V) = names(object$coef)
    V
}

#' Extract model residuals.
#'
#' @param object an object of class \code{sparse.sglmm}, typically the result of a call to \code{\link{sparse.sglmm}}.
#' @param type the type of residuals that should be returned. The alternatives are \dQuote{\code{deviance}} (default), \dQuote{\code{pearson}}, and \dQuote{\code{response}}.
#' @param \dots additional arguments.
#' @return A vector of residuals.
#' @seealso \code{\link{sparse.sglmm}}, \code{\link{residuals.glm}}
#' @method residuals sparse.sglmm
#' @export

residuals.sparse.sglmm = function(object, type = c("deviance", "pearson", "response"), ...)
{
    type = match.arg(type)
    if (type == "response")
        return(object$residuals)
    else if (type == "deviance")
    {
        mu.hat = object$fitted.values
        y = if (is.null(object$y))
                object$residuals + mu.hat
            else object$y
        wts = rep(1, length(y))
        d.res = sqrt(pmax((object$family$dev.resids)(y, mu.hat, wts), 0))
        return(ifelse(y > mu.hat, d.res, -d.res))
    }
    else # type == "pearson"
    {
        mu.hat = object$fitted.values
        se = sqrt(object$family$variance(mu.hat))
        return(object$residuals / se)
    }
}

#' Print a summary of a sparse SGLMM fit.
#'
#' @details This function displays (1) the call to \code{\link{sparse.sglmm}}, (2) the values of the hyperparameters and tuning parameters, (3) a table of estimates, (4) the DIC value for the fit, and (5) the number of posterior samples. Each row of the table of estimates shows an estimated regression coefficient, the HPD interval for the coefficient, and the Monte Carlo standard error.
#' @param object an object of class \code{sparse.sglmm}, typically the result of a call to \code{\link{sparse.sglmm}}.
#' @param alpha the significance level used to compute the HPD intervals. The default is 0.05.
#' @param digits the number of significant digits to display. The default is 4.
#' @param \dots additional arguments.
#' @seealso \code{\link{sparse.sglmm}}
#' @method summary sparse.sglmm
#' @export

summary.sparse.sglmm = function(object, alpha = 0.05, digits = 4, ...)
{
	if (! is.numeric(alpha) || length(alpha) > 1 || alpha <= 0 || alpha >= 1)
	    stop("'alpha' must be a number between 0 and 1.")
    cat("\nCall:\n\n")
    print(object$call)
    if (length(object$tune) > 0)
    {
        cat("\nTuning parameters:\n")
        tune.table = cbind(unlist(object$tune))
        colnames(tune.table) = ""
        print(tune.table, quote = FALSE)
    }
    if (length(object$hyper) > 0)
    {
        cat("\nHyperparameters:\n")
        hyper.table = cbind(unlist(object$hyper))
        colnames(hyper.table) = ""
        print(hyper.table, quote = FALSE)
    }
    p = length(object$coef)
    ci = matrix(, p, 2)
    mcse = object$beta.mcse
    for (j in 1:p)
        ci[j, ] = hpd(object$beta.sample[, j], alpha)
    coef.table = cbind(object$coef, ci, mcse)
    colnames(coef.table) = c("Estimate", "Lower", "Upper", "MCSE")
    rownames(coef.table) = names(object$coef)
    cat("\nCoefficients:\n\n")
    print(signif(coef.table, digits))
    cat("\nDIC:", signif(object$dic, digits), "\n\nNumber of iterations:", object$iter, "\n\n")
}

