#' Create adjacency matrix.
#'
#' @details This function builds the adjacency matrix for a \code{m} by \code{n} square lattice.
#'
#' @param m the number of rows in the lattice.
#' @param n the number of columns in the lattice. Defaults to \code{NULL}. If missing, the lattice is assumed to be \code{m} by \code{m}.
#'
#' @return A matrix \eqn{A} of 0s and 1s, where \eqn{(A_{ij})} is equal to 1 if and only if vertices \eqn{i} and \eqn{j} are adjacent.
#'
#' @export

adjacency.matrix = function(m, n = NULL)
{
    if (missing(n))
    {
        A = matrix(0, m^2, m^2)
        for (i in 1:m^2)
        {
            up = i - m
            down = i + m
            left = i - 1
            right = i + 1
            if (up > 0)
                A[i, up] = 1
            if (down <= m^2)
                A[i, down] = 1
            if (left %% m != 0)
                A[i, left] = 1
            if (i %% m != 0)
                A[i, right] = 1
        }
    }
    else
    {
        A = matrix(0, m * n, m * n)
        for (i in 1:(m * n))
        {
            up = i - n
            down = i + n
            left = i - 1
            right = i + 1
            if (up > 0)
                A[i, up] = 1
            if (down <= (m * n))
                A[i, down] = 1
            if (left %% n != 0)
                A[i, left] = 1
            if (i %% n != 0)
                A[i, right] = 1
        }
    }
    return(A)
}


build.M = function(A, d, rhoMax, epsilon)                        # compute M used to approx CAR variances.
{
    W = A / d
    cG = max(diag(W %*% W) / d)
    k =  ceiling((log(1 - rhoMax) + log(epsilon) - log(cG)) / log(rhoMax) - 1)
    eig = eigen(W)
    E = eig$vectors
    lam = eig$values
    Einv = solve(E)
    B = E * t(Einv)
    result = suppressWarnings(buildM$buildMcpp(B, k, lam))       # call C++ func from buildM Rcpp module.
    return(result)
}


neighbor.index = function(A)                                     # create index of adjacent neighbor pairs.
{
    index = c(NA, NA)                                            # initialize index.
    for (i in 1:nrow(A))
    {
        for (j in 1:ncol(A))
        {
            if (A[i, j] == 1 && i < j)
		index = rbind(index, c(i, j))                    # index (i, j) \in E with i < j.
	  }
    }
    index = index[-1, ]
    return(index)
}


#' Simulate areal data.
#'
#' @description \code{rcopCAR} simulates areal data from the copCAR model.
#'
#' @details This function randomly generates Poisson or Bernoulli areal data with adjacency matrix \eqn{A} from the copCAR model with the given spatial dependence parameter \eqn{\rho}, regression coefficients \eqn{\beta = (\beta_1, \dots, \beta_p)'}, and design matrix \eqn{X}. For more details on the copCAR model, see \code{\link{copCAR}}.
#'
#' @param rho the spatial dependence parameter \eqn{\rho} such that \eqn{\rho \in [0, 1)}.
#' @param beta the vector of regression coefficients  \eqn{\beta = (\beta_1, \dots, \beta_p)'}.
#' @param X the \eqn{n} by \eqn{p} design matrix \eqn{X}.
#' @param A the symmetric binary adjacency matrix for the underlying graph. \code{\link{adjacency.matrix}} generates an adjacency matrix for the \eqn{m} by \eqn{n} square lattice.
#' @param family the marginal distribution of the observations and link function to be used in the model. This can be a character string naming a family function, a family function or the result of a call to a family function. (See \code{\link{family}} for details of family functions.) Supported familes are \code{binomial} and \code{poisson}.
#'
#' @return A vector of length \eqn{n} distributed according to the copCAR model with the given design matrix and parameter values.
#'
#' @export
#'
#' @examples \dontrun{
#' require(lattice)
#'
#' # Use the 20 x 20 square lattice as the underlying graph.
#' m = 20
#' A = adjacency.matrix(m)
#'
#' # Set dependence parameter and regression coefficients.
#' rho = 0.8
#' beta = c(1, 1)
#'
#' # Create design matrix by assigning coordinates to each vertex
#' # such that the coordinates are restricted to the unit square.
#' x = rep(0:(m - 1) / (m - 1), times = m)
#' y = rep(0:(m - 1) / (m - 1), each = m)
#' X = cbind(x, y)
#'
#' # Draw Poisson data from copCAR model.
#' Z.pois = rcopCAR(rho, beta, X, A, family = poisson(link = "log"))
#'
#' # Create a level plot of the simulated data.
#' dev.new()
#' levelplot(Z ~ x * y, aspect = "iso")
#'
#' # Draw Bernoulli data from copCAR model.
#' Z.ber = rcopCAR(rho, beta, X, A, family = binomial(link = "logit"))
#'
#' # Create a level plot of the simulated data.
#' dev.new()
#' levelplot(Z2 ~ x * y, aspect = "iso")
#' }


rcopCAR = function(rho, beta, X, A, family)
{
    if (missing(family))
        stop("You must supply a family.")
    if (is.character(family))
        family = get(family, mode = "function", envir = parent.frame())
    if (is.function(family))
        family = family()
    if (is.null(family$family))
        stop("'family' not recognized.")
    if (! family$family %in% c("binomial", "poisson"))
        stop("'family' must be binomial or poisson.")
    if (missing(rho) || missing(beta))
        stop("You must supply the parameters rho and beta.")
    if (missing(A) || !is.matrix(A) || !isSymmetric(A) || !(A == 0 || A == 1))
        stop("You must supply a symmetric binary adjacency matrix.")
    if (rho <= 0 || rho >= 1)
        stop("'rho' must be between 0 and 1.")
    if (missing(X))
        stop("You must supply the design matrix X.")
    if (length(beta) != ncol(X))
        stop("The length of 'beta' must equal the number of columns of 'X'")
    Z = draw.copCAR(rho = rho,
                    beta = beta,
                    X = X,
                    A = A,
                    family = family)
    return(as.vector(Z))
}


draw.copCAR = function(rho, beta, X, A, family)
{
    d = rowSums(A)
    D = diag(d, length(d))
    Q = D - rho * A
    R = solve(Q)                                                 # find Q.inv.
    rootR = t(chol(R))
    Y = rootR %*% rnorm(length(d))                               # draw from MVN.
    v = diag(R)                                                  # CAR variances.
    U = pnorm(Y, 0, sqrt(v))                                     # prob integral transform (PIT).
    linkinv = family$linkinv
    eta = X %*% beta                                             # marginal linear predictors.
    mu = linkinv(eta)                                            # marginal means.
    if (family$family == "binomial")
    {
        Z = qbinom(U, 1, mu)                                     # inverse PIT.
    }
    if (family$family == "poisson")
    {
        Z = qpois(U, mu)                                         # inverse PIT.
    }
    return(Z)
}


copCAR.DT = function(Z, X, A, rhoMax, epsilon, offset, family)
{
    d = rowSums(A)
    M = build.M(A, d, rhoMax, epsilon)
    fit.glm = glm(Z ~ X - 1, offset = offset, family = family)   # get initial values for beta.
    D = diag(d, length(d))
    C = chol(as.spam(D - 0.5*A))                                 # Cholesky decomp of Q using initial value of rho.
    fit = optim(c(qnorm(0.5), fit.glm$coef),                     # optimize.
                lik.DT,
                Z = Z,
                X = X,
                A = A,
                D = D,
                d = d,
                C = C,
                M = M,
                offset = offset,
                family = family,
                method = "L-BFGS-B",
                lower = c(qnorm(1e-6), rep(-Inf, ncol(X))),
                upper = c(qnorm(0.999), rep(Inf, ncol(X))),
                hessian = TRUE)
    fit$par[1] = pnorm(fit$par[1])
    return(fit)
}


lik.DT = function(params, Z, X, A, D, d, C, M, offset, family)
{
    rho = pnorm(params[1])
    beta = params[-1]
    k = ncol(M) - 1
    v = as.vector(M %*% rho^(0:k) / d)                           # approx. CAR variances.
    Q = as.spam(D - rho * A)
    update = try(update.spam.chol.NgPeyton(C, Q), silent = TRUE) # update Cholesky decomp with new rho.
    linkinv = family$linkinv
    eta = offset + X %*% beta                                    # marginal linear predictor.
    mu = linkinv(eta)                                            # marginal means.
    U = ppois(Z, mu) - 0.5 * dpois(Z, mu)                        # apply DT.
    Y = qnorm(U, 0, sqrt(v))
    Y = ifelse(Y == Inf, 1e6, Y)                                 # deal with infinite values.
    qform =  t(Y) %*% (Q - diag(1/v, length(d))) %*% Y           # quadratic form for loglik.
    if (family$family == "poisson")
        loglik = - sum(log(diag(C))) - 0.5 * sum(log(v)) + 0.5 * qform - sum(dpois(Z, mu, log = TRUE))
    if (family$family == "binomial")
        loglik = - sum(log(diag(C))) - 0.5 * sum(log(v)) + 0.5 * qform - sum(dbinom(Z, 1, mu, log = TRUE))
    if (loglik == Inf)
        loglik = 1e6
    if (loglik == -Inf)
        loglik = -1e6
    return(loglik)
}


copCAR.CE = function(Z, X, A, rhoMax, epsilon, offset, family, m)
{
    d = rowSums(A)
    M = build.M(A, d, rhoMax, epsilon)
    U = matrix(runif(m * length(d)), length(d), m)               # draw standard uniforms.
    fit.glm = glm(Z ~ X - 1, offset = offset, family = family)   # get initial values for beta.
    D = diag(d, length(d))
    C = chol(as.spam(D - 0.5*A))                                 # Cholesky decomp of Q using initial value of rho.
    fit = optim(c(qnorm(0.5), fit.glm$coef),                     # optimize.
                lik.CE,
                Z = Z,
                X = X,
                A = A,
                D = D,
                d = d,
                C = C,
                M = M,
                U = U,
                offset = offset,
                family = family,
                method = "L-BFGS-B",
                lower = c(qnorm(1e-6), rep(-Inf, ncol(X))),
                upper = c(qnorm(0.999), rep(Inf, ncol(X))),
                hessian = TRUE)
    fit$par[1] = pnorm(fit$par[1])
    return(fit)
}


lik.CE = function(params, Z, X, A, D, d, C, M, U, offset, family)
{
    rho = pnorm(params[1])
    beta = params[-1]
    k = ncol(M) - 1
    v = as.vector(M %*% rho^(0:k) / d)                           # approx. CAR variances.
    Q = as.spam(D - rho * A)
    update = try(update.spam.chol.NgPeyton(C, Q), silent = TRUE) # update Cholesky decomp with new rho.
    linkinv = family$linkinv
    eta = offset + X %*% beta                                    # marginal linear predictor.
    mu = linkinv(eta)                                            # marginal means.
    m = ncol(U)
    w = numeric(m)                                               # initialize for w.
    for (j in 1:m)
    {
        z_ = Z - U[, j]                                          # continue Z.
        u_ = ppois.CE(z_, mu)
        y_ = qnorm(u_, 0, sqrt(v))
        qform = try(t(y_) %*% Q %*% y_, silent = TRUE)           # quadratic form for loglik.
        if (class(qform) == "try-error")                         # deal with errors in calc of quadratic form.
            return(10000)
        w[j] = -0.5 * (qform - sum(y_^2 / v))
    }
    s = max(w)                                                   # choose s for log-sum-exp trick.
    loglik = -log(sum(exp(w - s))) - s - sum(dpois(Z, mu, log = TRUE)) - sum(log(diag(C))) - 0.5 * sum(log(v))
    if (loglik == Inf)
       loglik = 1e6
    if (loglik == -Inf)
       loglik = -1e6
    return(loglik)
}


dpois.CE = function(x, mu)
{
    dpois(floor(x + 1), mu)                                      # pdf for continued Z.
}


ppois.CE = function(q, mu)
{
    ppois(floor(q), mu) + (q - floor(q)) * dpois.CE(q, mu)       # cdf for continued Z.
}


compute.mcse = function(Z, X, A, rhoMax, epsilon, offset, family, m, mcse.iter, verbose)
{
    rho.hat = numeric(mcse.iter)                                 # initialze vector for rho.
    for (k in 1:mcse.iter)
    {
        if (verbose == TRUE)
        {
            cat("Progress => ", round(k / mcse.iter * 100, 0), "%\n", sep = "")
            flush.console()
        }
        fit = copCAR.CE(Z = Z,                                    # fit the model mcse.iter times.
                        X = X,
                        A = A,
                        rhoMax = rhoMax,
                        epsilon = epsilon,
                        offset = offset,
                        family = family,
                        m = m)
        rho.hat[k] = fit$par[1]
    }
    mcse = sum((rho.hat - mean(rho.hat))^2)
    mcse = sqrt(mcse / (mcse.iter - 1))                          # compute mcse for rho.
    cv = mcse / mean(rho.hat)                                    # compute cv for rho.
    result = list()
    result$mcse = mcse                                           # return mcse.
    result$cv = cv                                               # return cv.
    return(result)
}


copCAR.CML = function(Z, X, A, offset, family)
{
    d = rowSums(A)
    D = diag(d, length(d))
    index = neighbor.index(A)                                    # get index of adjacent pairs.
    fit.glm = glm(Z ~ X - 1, offset = offset, family = family)   # get initial values for beta.
    fit = optim(c(qnorm(0.5), fit.glm$coef),                     # optimize.
                lik.CML,
                Z = Z,
                X = X,
                A = A,
                D = D,
                d = d,
                offset = offset,
                family = family,
                index = index,
                method = "L-BFGS-B",
                lower = c(qnorm(1e-6), rep(-Inf, ncol(X))),
                upper = c(qnorm(0.999), rep(Inf, ncol(X))),
                hessian = TRUE)
    fit$par[1] = pnorm(fit$par[1])
    return(fit)
}


lik.CML = function(params, Z, X, A, D, d, offset, family, index)
{
    rho = pnorm(params[1])
    beta = params[-1]
    linkinv = family$linkinv
    eta = offset + X %*% beta                                    # marginal linear predictors.
    mu = linkinv(eta)                                            # marginal means.
    Q = D - rho * A
    R = inverse$armaInv(Q)                                       # find Q.inv: call C++ func from inverse Rcpp module.
    v = diag(R)                                                  # CAR variances.
    if (family$family == "binomial")                             # values for bivariate copula.
    {
        Y.0 = qnorm(pbinom(Z, 1, mu), 0, sqrt(v)) / sqrt(v)
        Y.1 = qnorm(pbinom(Z - 1, 1, mu), 0, sqrt(v)) / sqrt(v)
    }
    if (family$family == "poisson")                              # values for bivariate copula.
    {
        Y.0 = qnorm(ppois(Z, mu), 0, sqrt(v)) / sqrt(v)
        Y.1 = qnorm(ppois(Z - 1, mu), 0, sqrt(v)) / sqrt(v)
    }
    Y.0 = ifelse(Y.0 == -Inf, -1e6, Y.0)                         # deal with nonfinite values.
    Y.1 = ifelse(Y.1 == -Inf, -1e6, Y.1)
    Y.0 = ifelse(Y.0 == Inf, 1e6, Y.0)
    Y.1 = ifelse(Y.1 == Inf, 1e6, Y.1)
    loglik = 0                                                   # initialize loglik.
    for (h in 1:nrow(index))                                     # index over adjacent pairs.
    {
        i = index[h, 1]
	  j = index[h, 2]
        corr = R[i, j] / sqrt((v[i] * v[j]))                     # correlation for bivariate copula.
        K = matrix(c(Y.0[i], Y.0[j],
                     Y.1[i], Y.0[j],
                     Y.0[i], Y.1[j],
                     Y.1[i], Y.1[j]),
                   nrow = 4,
                   byrow = TRUE)
        temp = pbivnormf(K, corr)                                # bivariate copula for each pair of obs.
        if (sum(temp[-(2:3)] - temp[2:3]) < 0) {
            temp = 0
        } else {
            temp = - log(sum(temp[-(2:3)] - temp[2:3])) }
	  if (temp == Inf)
            temp = 1e6
	  if (temp == -Inf)
	    temp = -1e6
        loglik = loglik + temp                                   # sum log likelihood.
    }
    if (loglik == Inf)
       loglik = 1e6
    if (loglik == -Inf)
       loglik = -1e6
    return(loglik)
}


pbivnormf = function (K, rho = 0)                                # spec of pbivnorm wrapper.
{
    correl = rep(rho, nrow(K))
    lower = as.double(c(0, 0))
    infin = as.integer(c(0, 0))
    uppera = as.double(K[, 1])
    upperb = as.double(K[, 2])
    lt = as.integer(nrow(K))
    prob = double(lt)
    correl = as.double(correl)
    ans = .Fortran("PBIVNORM", prob, lower, uppera, upperb,
        infin, correl, lt, PACKAGE = "copCAR")[[1]]
    return(ans)
}


CI.asymptotic = function(fit, conf.level)
{
    crit.val = qnorm((1 + conf.level) / 2)
    params = fit$par
    theta = c(qnorm(params[1]), params[-1])                      # reparameterization.
    npar = length(theta)                                         # number of parameters
    I.hat.inv = solve(fit$hessian)                               # inverse obs Fisher info
    CI = matrix(NA, npar, 2)
    for (i in 1:npar)
        CI[i, ] = theta[i] + c(-1, 1) * crit.val * sqrt(I.hat.inv[i, i])
    CI[1, ] = pnorm(CI[1, ])                                     # undo reparameterization.
    colnames(CI) = c("lower CI", "upper CI")
    rownames(CI) = c("rho", paste("beta", 1:(npar-1), sep = ""))
    result = list()
    result$CI = CI                                               # return CI.
    result$varcov = I.hat.inv                                    # return inverse obs Fisher info.
    return(result)
}


CI.bootstrap = function(fit, X, A, b, conf.level, offset, family, type, verbose, rhoMax, epsilon)
{
    crit.val = qnorm((1 + conf.level) / 2)
    params = fit$par
    theta = c(qnorm(params[1]), params[-1])                      # reparameterization.
    npar = length(params)                                        # number of parameters
    d = rowSums(A)
    D = diag(d, length(d))
    Q = D - params[1] * A
    R = solve(Q)                                                 # calculate Q.inv.
    rootR = t(chol(R))
    v = diag(R)                                                  # CAR variances.
    J.hat = matrix(0, nrow = npar, ncol = npar)                  # initialize J.hat matrix.
    linkinv = family$linkinv
    eta = offset + X %*% params[-1]                              # marginal linear predictors.
    mu = linkinv(eta)                                            # marginal means.
    if (type == "DT")
    {
        M = build.M(A, d, rhoMax, epsilon)                       # calc M, C for gradient
        C = chol(as.spam(Q))
    }
    if (type == "CML")
        index = neighbor.index(A)                                # get index of adjacent pairs.
    for (j in 1:b)                                               # parametric bootstrap: b iterations.
    {
        if (verbose == TRUE)
        {
            cat("Progress => ", round(j / b * 100, 0), "%\n", sep = "")
            flush.console()
        }
        Z_ = boot.sim.helper(rootR = rootR,                      # draw Z_.
                             d = d,
                             v = v,
                             family = family,
                             mu = mu)
        if (type == "CML")
        {
            gr = - grad(lik.CML,                                 # calculate gradient.
                        theta,
                        Z = Z_,
                        X = X,
                        A = A,
                        D = D,
                        d = d,
                        offset = offset,
                        family = family,
                        index = index)
        }
        if (type == "DT")
        {
            gr = - grad(lik.DT,                                  # calculate gradient.
                        theta,
                        Z = Z_,
                        X = X,
                        A = A,
                        D = D,
                        d = d,
                        C = C,
                        M = M,
                        offset = offset,
                        family = family)
         }
         J.hat = J.hat + gr %o% gr / b                           # get the tempeh...
    }
    I.hat.inv = solve(fit$hessian)                               # then get the bread...
    G.hat.inv = I.hat.inv %*% J.hat %*% I.hat.inv                # and make the sandwich.
    CI = matrix(NA, npar, 2)
    for (i in 1:npar)
	CI[i, ] = theta[i] + c(-1, 1) * crit.val * sqrt(G.hat.inv[i, i])
    CI[1, ] = pnorm(CI[1, ])                                     # undo reparameterization.
    colnames(CI) = c("lower CI", "upper CI")
    rownames(CI) = c("rho", paste("beta", 1:(npar-1), sep = ""))
    result = list()
    result$CI = CI                                               # return CI.
    result$varcov = G.hat.inv                                    # return sandwich.
    return(result)
}


boot.sim.helper = function(rootR, d, v, family, mu)
{
    Y = rootR %*% rnorm(length(d))                               # draw from MVN.
    U = pnorm(Y, 0, sqrt(v))                                     # prob integral trans (PIT).
    if (family$family == "binomial")
    {
        Z = qbinom(U, 1, mu)                                     # inverse PIT.
    }
    if (family$family == "poisson")
    {
        Z = qpois(U, mu)                                         # inverse PIT.
    }
    return(Z)
}


#' Fit copCAR model to discrete areal data.
#'
#' @description Fit the copCAR model to areal data consisting of Poisson or Bernoulli marginal observations.
#'
#' @details This function performs frequentist inference for the copCAR model of Hughes (2014), a copula-based areal regression model that uses the conditional autoregression (CAR) from the spatial generalized linear mixed model (Besag, 1974). Specifically, copCAR uses the CAR copula, a Caussian copula based on the proper CAR. The CAR copula is specified as \deqn{C_{Q^{-1}}(u) = \Phi_{Q^{-1}}(\Phi_{\sigma_1}^{-1}(u_1), \dots, \Phi_{\sigma_n}^{-1}(u_n)),} where \eqn{\Phi_{\sigma_i}} denotes the cdf of the normal distribution with mean zero and variance \eqn{\sigma_i^2}, \eqn{Q = D - \rho A} such that \eqn{\tau Q} is the precision matrix of the proper CAR, \eqn{A} is the adjacency matrix for the underlying graph, \eqn{D = diag(d_1, \dots, d_n)} where \eqn{d_i} is the degree of vertex \eqn{i} of the underlying graph, and \eqn{u = (u_1, \dots, u_n)'} is a realization of the copula such that \eqn{z_i = F_i^{-1}(u_i)} for the marginal observation \eqn{z_i} having desired marginal distribution function \eqn{F_i}. For Bernoulli marginals, the expectation is \eqn{(1 + \exp(-x_i'\beta))^{-1}}; for Poisson marginals, the expectation is \eqn{\exp(x_i'\beta)}, where \eqn{\beta = (\beta_1, \dots, \beta_p)'} is the regression coefficient. Note that the CAR variances \eqn{(\sigma_1^2, \dots, \sigma_n^2)' = vecdiag(Q^{-1})} are not free parameters but are determined by the spatial dependence parameter \eqn{\rho}.
#' \cr
#' \cr
#' The spatial dependence parameter \eqn{\rho \in [0, 1)} and regression coefficient \eqn{\beta = (\beta_1, \dots, \beta_p)' \in R^p} can be estimated using the continous extension (CE) (Madsen, 2009), distribtional transform (DT) (Kazianka and Pilz, 2010), or composite marginal likelihood (CML) (Varin, 2008).
#' \cr
#' \cr
#' The CE approach optimizes an approximate maximum likelihood by sampling \eqn{m} independent standard uniform vectors of length \eqn{n} used to transform the discrete observations into continous random variables via convolution (Denuit and Lambert, 2005). The size of \eqn{m} can be choosen by computing Monte Carlo standard errors (Flegal et al., 2008). If the Monte Carlo standard error of the estimate for \eqn{\rho} is small relative to the sample mean, that is, if the estimated coefficient of variation is sufficiently small, the current value of \eqn{m} is sufficiently large. The CE is exact up to Monte Carlo standard error, but is computationally intensive and not suitable for Bernoulli marginals. If requested, asymptotic confidence intervals for the parameters are computed using the observed inverse Fisher information.
#' \cr
#' \cr
#' The DT stochastically "smoothes" the jumps of the discrete distribution function, an approach that goes at least as far back as Ferguson (1967). The DT-based approximation performs well for Poisson marginals. Since the log-likelihood is misspecified, the asympototic covariance matrix is the Godambe information matrix (Godambe, 1960). This is estimated using a parametric bootstrap for the variance of the score when computing confidence intervals for the parameters.
#' \cr
#' \cr
#' The CML approach specifies the likelihood as a product of pairwise likelihoods of adjacent observations, and performs well for both Poisson and Bernoulli data. Similar to the DT, the log-likelihood is misspecified, so the confidence intervals for the parameters are computed via a parametric bootstrap.
#' \cr
#' \cr
#' In the CE and DT approaches, the CAR variances are approximated by \eqn{(\tilde{\sigma}_1^2, \dots, \tilde{\sigma}_n^2)'} such that \eqn{(\sigma_i^2 - \tilde{\sigma}_i^2) < \epsilon} for every \eqn{i = 1, \dots, n} for a specified tolerance \eqn{\epsilon > 0} and every \eqn{\rho \in [0, \rho^{\max})}.
#' \cr
#' \cr
#' @param formula an object of class "\code{\link{formula}}" (or one that can be coerced to that class): a symbolic description of the model to be fitted. The details of the model specification are given under "Details".
#' @param A the symmetric binary adjacency matrix for the underlying graph. \code{\link{adjacency.matrix}} generates an adjacency matrix for the \eqn{m} by \eqn{n} square lattice.
#' @param family the marginal distribution of the observations at the areal units and link function to be used in the model. This can be a character string naming a family function, a family function or the result of a call to a family function. (See \code{\link{family}} for details of family functions.) Supported families are \code{binomial} and \code{poisson}.
#' @param method the method for inference. \code{copCAR} supports the continous extension ("\code{CE}"), distributional transform ("\code{DT}"), and composite marginal likelihood ("\code{CML}").
#' @param conf.int the method for computing confidence intervals. \code{"asympototic"} is appropriate for the continous extension. \code{"bootstrap"} performs a parametric boostrap appropriate for the distributional transform and composite marginal likelihood.
#' @param data an optional data frame, list or environment (or object coercible by \code{\link{as.data.frame}} to a data frame) containing the variables in the model. If not found in data, the variables are taken from \code{environment(formula)}, typically the environment from which \code{copCAR} is called.
#' @param offset this can be used to specify an \emph{a priori} known component to be included in the linear predictor during fitting. This should be \code{NULL} or a numeric vector of length equal to the number of observations. One or more \code{\link{offset}} terms can be included in the formula instead or as well, and if more than one is specified their sum is used. See \code{\link{model.offset}}.
#' @param control a list of parameters for controlling the fitting process.
#' \describe{
#'        \item{\code{conf.level}}{the value \eqn{1 - \alpha} used for computing confidence intervals. Defaults to \code{0.95}.}
#'        \item{\code{boot.iter}}{the size of the parametric bootstrap sample. This applies when \code{conf.int = "bootstrap"}. Defaults to \code{500}.}
#'        \item{\code{m}}{the number of independent standard uniforms used to approximate the expected likelhood when \code{method = "CE"}. This applies when \code{method = "CE"}. Defaults to \code{1000}.}
#'        \item{\code{mcse}}{logical. Should the Monte Carlo standard error for \eqn{\rho} be computed? Use only when \code{method = "CE"}. Defaults to FALSE. The Monte Carlo standard error is calculated using a sample size of \code{mcse.iter}.}
#'        \item{\code{mcse.iter}}{the Monte Carlo standard error sample size for \eqn{\rho} when \code{method = "CE"} and \code{mcse = TRUE}. Defaults to \code{100}.}
#'        \item{\code{rhoMax}}{the value \eqn{\rho^{\max}}, which is the maximum value of \eqn{\rho} used to approximate the CAR variances when \code{method = "CE"} or \code{method = "DT"}. If missing, assumed to be \code{0.999}.}
#'        \item{\code{epsilon}}{the tolerance \eqn{\epsilon > 0} used to approximate the CAR variances when \code{method = "CE"} or \code{method = "DT"}. If missing, assumed to be \code{0.001}.}
#'        \item{\code{verbose}}{a logical value indicating whether to print bootstrap or mcse progress to the screen. Defaults to \code{FALSE}.}
#' }
#'
#' @return \code{copCAR} returns an object of S3 class \code{"copCAR"}, which is a list containing the following components:
#'         \item{coefficients}{the point estimate of \eqn{(\rho, \beta')'}.}
#'         \item{conf.int}{(if \code{conf.int} is not \code{"none"}) the confidence intervals for \eqn{(\rho, \beta')'}.}
#'         \item{conf.type}{ the type of confidence interval specified.}
#'         \item{conf.level}{(if \code{conf.int} is not \code{"none"}) the confidence level, \eqn{1 - \alpha}.}
#'         \item{mcse}{(if \code{method = "CE"} and \code{mcse = TRUE}) the Monte Carlo standard error of \eqn{\rho}.}
#'         \item{mcse.iter}{(if \code{method = "CE"} and \code{mcse = TRUE}) the Monte Carlo standard error sample size.}
#'         \item{mcse.cv}{(if \code{method = "CE"} and \code{mcse = TRUE}) the estimated coefficient of variation of \eqn{\rho}.}
#'         \item{I.inv}{(if \code{conf.int = "asymptotic"}) the estimated inverse observed Fisher information matrix (hessian) for \eqn{(\Phi^{-1}(\rho), \beta')'}.}
#'         \item{G.inv}{(if \code{conf.int = "bootstrap"}) the estimated inverse Godambe information matrix (sandwich) for \eqn{(\Phi^{-1}(\rho), \beta')'}.}
#'         \item{se}{(if \code{conf.int = "bootstrap"} or \code{conf.int = "asymptotic"}) the estimated standard errors for \eqn{(\Phi^{-1}(\rho), \beta')'}.}
#'         \item{boot.iter}{(if \code{conf.int = "bootstrap"}) the number of parametric bootstrap samples.}
#'         \item{Z}{the response vector used.}
#'         \item{X}{the design matrix.}
#'         \item{model}{the model frame.}
#'         \item{npar}{the number of model parameters.}
#'         \item{marginal.linear.predictors}{linear predictors for the margins.}
#'         \item{marginal.fitted.values}{fitted values for the margins.}
#'         \item{call}{the matched call.}
#'         \item{formula}{the formula supplied.}
#'         \item{method}{the method used for inference.}
#'         \item{convergence}{the integer code returned by \code{\link{optim}} subsequent to optimizing the log-likelihood.}
#'         \item{message}{a character string to go along with \code{convergence}.}
#'         \item{terms}{the \code{terms} object used.}
#'         \item{data}{the \code{data} argument.}
#'         \item{xlevels}{(where relevant) a record of the levels of the factors used in fitting.}
#'         \item{control}{a list containing the names and values of the control parameters.}
#'         \item{value}{the value of the negative log-likelihood at its minimum.}
#'
#' @references
#' Besag, J. (1974) Spatial interaction and the statistical analysis of lattice systems.  \emph{Journal of the Royal Statistical Soceity, Series B, Methodological}, 36(2), 192--236.
#' @references
#' Denuit, M. and Lambert, P. (2005) Constraints on concordance measures in bivariate discrete data. \emph{Journal of Multivariate Analysis}, 93, 40--57.
#' @references
#' Ferguson, T. (1967) \emph{Mathematical statistics: a decision theoretic approach}, New York: Academic Press.
#' @references
#' Flegal, J., Haran, M., and Jones, G. (2008) Markov Chain Monte Carlo: can we trust the third significant figure? \emph{Statistical Science}, 23(2), 250--260.
#' @references
#' Godambe, V. (1960) An optimum property of regular maximum likelihood estimation. \emph{The Annals of Mathmatical Statistics}, 31(4), 1208--1211.
#' @references
#' Kazianka, H. and Pilz, J. (2010) Copula-based geostatistical modeling of continuous and discrete data including covariates. \emph{Stochastic Environmental Research and Risk Assessment}, 24(5), 661--673.
#' @references
#' Madsen, L. (2009) Maximum likelihood estimation of regression parameters with spatially dependent discrete data. \emph{Journal of Agricultural, Biological, and Environmental Statistics}, 14(4), 375--391.
#' @references
#' Varin, C. (2008) On composite marginal likelihoods. \emph{Advances in Statistical Analysis}, 92(1), 1--28.
#'
#' @export
#'
#' @examples \dontrun{
#' # Simulate data and fit copCAR model.
#'
#' # Use the 20 x 20 square lattice as the underlying graph.
#' m = 20
#' A = adjacency.matrix(m)
#'
#' # Set dependence parameter and regression coefficients.
#' rho = 0.8
#' beta = c(1, 1)
#'
#' # Create design matrix by assigning coordinates to each vertex
#' # such that the coordinates are restricted to the unit square.
#' x = rep(0:(m - 1) / (m - 1), times = m)
#' y = rep(0:(m - 1) / (m - 1), each = m)
#' X = cbind(x, y)
#'
#' # Draw Poisson data from copCAR model.
#' Z = rcopCAR(rho, beta, X, A, family = poisson(link = "log"))
#'
#' # Fit the copCAR model using the continous extension and
#' # compute 95% (default) aysmptotic CI for rho and beta.
#' fit.CE = copCAR(Z ~ X - 1, A, family = poisson, method = "CE", conf.int = "asymptotic")
#' summary(fit.CE)
#'
#' # Fit the copCAR model using the distributional transform and
#' # compute 90% CI for rho and beta using 100 bootstrap iterations.
#' fit.DT = copCAR(Z ~ X - 1, A, family = poisson, method = "DT", conf.int = "bootstrap",
#'                 control = list(conf.level = 0.90, boot.iter = 100))
#' summary(fit.DT)
#'
#' # Fit the copCAR model using composite marginal likelihood and
#' # do not compute a CI for rho and beta.
#' fit.CML = copCAR(Z ~ X - 1, A, family = poisson, method = "CML", conf.int = "none")
#' summary(fit.CML)
#' }


copCAR = function(formula,
                  A,
                  family,
                  method = c("CML", "DT", "CE"),
                  conf.int = c("none", "bootstrap", "asymptotic"),
                  data,
                  offset = NULL,
                  control = list())
{
    call = match.call()
    if (missing(formula))
        stop("You must supply a formula.")
    if (is.character(family))
        family = get(family, mode = "function", envir = parent.frame())
    if (is.function(family))
        family = family()
    if (is.null(family$family))
        stop("'family' not recognized.")
    if (!family$family %in% c("binomial", "poisson"))
        stop("'family' must be binomial or poisson.")
    if (missing(data))
        data = environment(formula)
    if (missing(method))
        stop("You must supply a method.")
    if (!method %in% c("CML", "DT", "CE"))
        stop("Method not recognized.")
    if (missing(conf.int))
        conf.int = "none"
    if (!conf.int %in% c("none", "bootstrap", "asymptotic"))
        stop("Method not recognized.")
    if (missing(A) || !is.matrix(A) || !isSymmetric(A) || !(A == 0 || A == 1))
        stop("You must supply a symmetric binary adjacency matrix.")
    if (method != "CE" && conf.int == "asymptotic")
        stop("Confidence intervals for the DT and CML must be bootstrapped.")
    if (method == "CE" && conf.int == "bootstrap")
        stop("Bootsrap option for confidence intervals not available for the CE.")
    if (family$family == "binomial" && method == "CE")
        stop("The CE method cannot be used for Bernoulli observations.")
    if (family$family == "binomial" && method == "DT" && sqrt(nrow(A)) < 30 && sqrt(ncol(A)) < 30 )
        stop("The DT method performs poorly for Bernoulli observations with an underlying graph of less than 30 by 30.")
    if (!is.list(control))
        stop("'control' must be a list.")
    con = list(conf.level = 0.95, boot.iter = 500, m = 1000, mcse = FALSE, mcse.iter = 100, rhoMax = 0.999, epsilon = 0.001, verbose = FALSE)
    nmsC = names(con)
    namc = names(control)
    con[nmsC = names(control)] = control                         # set all missing control args to default.
    noNms = namc[!namc %in% nmsC]
    if (length(noNms) > 0)
        warning("Ignoring unknown names in control: ", paste(noNms, collapse = ", "))
    if (!is.logical(con$mcse))
        stop("'control$mcse' must be a logical value.")
    if (con$boot.iter <= 0 || con$boot.iter != trunc(con$boot.iter))
        stop("boot.iter must be a positive integer")
    if (con$mcse.iter <= 0 || con$mcse.iter != trunc(con$mcse.iter))
        stop("boot.iter must be a positive integer")
    if (con$m <= 0 || con$m != trunc(con$m))
        stop("boot.iter must be a positive integer")
    if (con$conf.level <= 0 || con$conf.level >= 1)
        stop("conf.level must be between 0 and 1")
    if (con$rhoMax <= 0 || con$rhoMax >= 1)
        stop("'rhoMax' must be between 0 and 1.")
    if (con$epsilon <= 0)
        stop("'epsilon' must be a positive number")
    if (!is.logical(con$verbose))
        stop("'verbose' must be a logical value.")
    mf = match.call(expand.dots = FALSE)
    m = match(c("formula", "data", "offset"), names(mf), 0L)
    mf = mf[c(1L, m)]
    mf[[1L]] = quote(stats::model.frame)
    mf = eval(mf, parent.frame())
    mt = attr(mf, "terms")
    Z = model.response(mf, "numeric")
    X = model.matrix(mt, mf)
    offset = as.vector(model.offset(mf))
    if (!is.null(offset) && length(offset) != length(Z))
        stop("number of offsets should equal number of observations.")
    if (is.null(offset))
        offset = rep(0, length(Z))
    if (sum(c(nrow(X), nrow(A)) != length(Z)) > 0)
        stop("The supplied response vector/design matrix/adjacency matrix are not conformable.")
    method = match.arg(method)
    if (method == "CE")                                          # fit if CE.
        fit = copCAR.CE(Z = Z,
                        X = X,
                        A = A,
                        rhoMax = con$rhoMax,
                        epsilon = con$epsilon,
                        offset = offset,
                        family = family,
                        m = con$m)
    if (method == "DT")                                          # fit if DT.
        fit = copCAR.DT(Z = Z,
                        X = X,
                        A = A,
                        rhoMax = con$rhoMax,
                        epsilon = con$epsilon,
                        offset = offset,
                        family = family)
    if (method == "CML")                                         # fit if CML.
        fit = copCAR.CML(Z = Z,
                         X = X,
                         A = A,
                         offset = offset,
                         family = family)
    if (is.null(fit$convergence) || fit$convergence != 0)        # check if optim converged.
        stop("Maximum likelihood estimation failed.")
    npar = length(fit$par)                                       # number of model parameters.
    linkinv = family$linkinv
    eta = X %*% fit$par[-1]                                      # marginal linear predictors.
    mu = linkinv(eta)                                            # marginal means.
    result = list()
    class(result) = "copCAR"                                     # create class copCAR for returned object.
    if (conf.int == "asymptotic")                                # compute CI if asymptotic.
    {
        CI = CI.asymptotic(fit, con$conf.level)
        result$I.inv = CI$varcov
        result$se = sqrt(diag(result$I.inv))                     # SE of theta.
        result$conf.level = con$conf.level
        result$conf.int = CI$CI
        names(result$se) = c("Phi.inv(rho)", paste("beta", 1:(npar-1), sep = ""))
        rownames(result$I.inv) = colnames(result$I.inv) = c("Phi.inv(rho)", paste("beta", 1:(npar-1), sep = ""))
    }
    if (conf.int == "bootstrap")                                 # compute CI if bootstrap.
    {
        CI = CI.bootstrap(fit,
                          X = X,
                          A = A,
                          b = con$boot.iter,
                          conf.level = con$conf.level,
                          offset = offset,
                          family = family,
                          type = method,
                          verbose = con$verbose,
                          rhoMax = con$rhoMax,
                          epsilon = con$epsilon)
        result$boot.iter = con$boot.iter
        result$G.inv = CI$varcov
        result$se = sqrt(diag(result$G.inv))                     # SE of theta.
        result$conf.level = con$conf.level
        result$conf.int = CI$CI
        names(result$se) = c("Phi.inv(rho)", paste("beta", 1:(npar-1), sep = ""))
        rownames(result$G.inv) = colnames(result$G.inv) = c("Phi.inv(rho)", paste("beta", 1:(npar-1), sep = ""))
    }
    if (method == "CE" && con$mcse == TRUE)                      # msce for rho in CE.
    {
        mcse = compute.mcse(Z = Z,
                            X = X,
                            A = A,
                            rhoMax = con$rhoMax,
                            epsilon = con$epsilon,
                            family = family,
                            offset = offset,
                            m = con$m,
                            mcse.iter = con$mcse.iter,
                            verbose = con$verbose)
        result$mcse = mcse$mcse
        result$mcse.iter = con$mcse.iter
        result$mcse.cv = mcse$cv
    }
    result$marginal.linear.predictors = eta
    result$marginal.fitted.values = mu
    result$npar = npar
    result$coefficients = fit$par
    result$conf.type = conf.int
    result$Z = Z
    result$X = X
    result$convergence = fit$convergence
    result$message = fit$message
    result$value = fit$value
    result$data = data
    result$method = method
    result$xlevels = .getXlevels(mt, mf)
    result$call = call
    result$terms = mt
    result$model = mf
    result$formula = formula
    con.out = list()                                             # control list to be returned.
    if (con$mcse == TRUE && method == "CE")
        con.out$mcse.iter = con$mcse.iter
    if (method == "CE")
        con.out$m = con$m
    if (conf.int != "none")
        con.out$conf.level = con$conf.level
    if (conf.int == "bootstrap")
        con.out$boot.iter = con$boot.iter
    result$control = con.out
    names(result$coefficients) = c("rho", paste("beta", 1:(npar - 1), sep = ""))
    return(result)
}


#' Print a summary of a copCAR model fit.
#'
#' @details This function displays (1) the call to \code{\link{copCAR}}, (2) the values of the control parameters, (3) a table of estimates, and (4) confidence intervals.
#'
#' Each row of the table of estimates shows a parameter estimate, the confidence interval for the parameter, and, where applicable, the Monte Carlo standard error.
#'
#' @param object an object of class \code{copCAR}, the result of a call to \code{\link{copCAR}}.
#' @param \dots additional arguments.
#'
#' @seealso \code{\link{copCAR}}
#'
#' @method summary copCAR
#'
#' @export

summary.copCAR = function(object, ...)
{
    cat("\nCall:\n\n")
    print(object$call)

    cat("\nControl parameters:\n")
    if (length(object$control) > 0)
    {
        control.table = cbind(unlist(object$control))
        colnames(control.table) = ""
        print(control.table, quote = FALSE)
    } else
    {
        print("None specified." )
    }
    npar = object$npar
    if (object$conf.type == "none")
        object$conf.int = matrix(rep(NA, 2 * npar), ncol = 2)
    coef.table = cbind(object$coefficients, object$conf.int)
    colnames(coef.table ) = c("Estimate", "Lower", "Upper")
    rownames(coef.table ) = c("rho", paste("beta", 1:(npar-1), sep = ""))
    if (object$method == "CE")
    {
        if (!is.null(object$mcse)) {
            mcse  = c(object$mcse, rep(NA, npar-1))
        } else {
            mcse = rep(NA, npar) }
        coef.table = cbind(coef.table, mcse)
    }
    cat("\nCoefficients:\n\n")
    print(signif(coef.table, digits = 4))
    cat("\nConvergence:\n\n")
    if (object$convergence == 0)
        print(paste("Optimization converged at loglik = -", signif(object$value, digits = 4), sep = ""))
    if (object$convergence != 0)
        print("Optimization did not converge.")
}


#' Return the covariance matrix of the parameters of a \code{copCAR} model object.
#'
#' @param object a fitted \code{copCAR} model object.
#' @param \dots additional arguments.
#'
#' @return An estimate of the covariance matrix of \eqn{(\Phi^{-1}(\rho), \beta')'}.
#'
#' @seealso \code{\link{copCAR}}
#'
#' @method vcov copCAR
#'
#' @export

vcov.copCAR = function(object, ...)
{
    if (object$conf.type == "none")
        stop("To compute the covariance matrix, the copCAR object needs to have confidence intervals calculated.")
    if (!is.null(object$I.inv))
        V = object$I.inv
    if (!is.null(object$G.inv))
        V = object$G.inv
    return(V)
}
