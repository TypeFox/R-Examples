
#' Compute Monte Carlo standard errors for expectations.
#'
#' @param x a vector of values from a Markov chain.
#' @param size the batch size. The default value is \dQuote{\code{sqroot}}, which uses the square root of the sample size. \dQuote{\code{cuberoot}} will cause the function to use the cube root of the sample size. A numeric value may be provided if neither \dQuote{\code{sqroot}} nor \dQuote{\code{cuberoot}} is satisfactory.
#' @param g a function such that \eqn{E(g(x))} is the quantity of interest. The default is \code{NULL}, which causes the identity function to be used.
#' @param method the method used to compute the standard error. This is one of \dQuote{\code{bm}} (batch means, the default), \dQuote{\code{obm}} (overlapping batch means), \dQuote{\code{tukey}} (spectral variance method with a Tukey-Hanning window), or \dQuote{\code{bartlett}} (spectral variance method with a Bartlett window).
#' @param warn a logical value indicating whether the function should issue a warning if the sample size is too small (less than 1,000).
#' @return \code{mcse} returns a list with two elements:
#'         \item{est}{an estimate of \eqn{E(g(x))}.}
#'         \item{se}{the Monte Carlo standard error.}
#' @references
#' Flegal, J. M. (2012) Applicability of subsampling bootstrap methods in Markov chain Monte Carlo. In Wozniakowski, H. and Plaskota, L., editors, \emph{Monte Carlo and Quasi-Monte Carlo Methods 2010} (to appear). Springer-Verlag.
#'
#' Flegal, J. M. and Jones, G. L. (2010) Batch means and spectral variance estimators in Markov chain Monte Carlo. \emph{The Annals of Statistics}, \bold{38}, 1034--1070.
#'
#' Flegal, J. M. and Jones, G. L. (2011) Implementing Markov chain Monte Carlo: Estimating with confidence. In Brooks, S., Gelman, A., Jones, G. L., and Meng, X., editors, \emph{Handbook of Markov Chain Monte Carlo}, pages 175--197. Chapman & Hall/CRC Press.
#'
#' Flegal, J. M., Jones, G. L., and Neath, R. (2012) Markov chain Monte Carlo estimation of quantiles. \emph{University of California, Riverside, Technical Report}.
#'
#' Jones, G. L., Haran, M., Caffo, B. S. and Neath, R. (2006) Fixed-width output analysis for Markov chain Monte Carlo. \emph{Journal of the American Statistical Association}, \bold{101}, 1537--1547.
#' @seealso
#' \code{\link{mcse.mat}}, which applies \code{mcse} to each column of a matrix or data frame.
#'
#' \code{\link{mcse.q}} and \code{\link{mcse.q.mat}}, which compute standard errors for quantiles.
#' @examples
#'
#' # Create 10,000 iterations of an AR(1) Markov chain with rho = 0.9.
#'
#' n = 10000
#' x = double(n)
#' x[1] = 2
#' for (i in 1:(n - 1))
#'     x[i + 1] = 0.9 * x[i] + rnorm(1)
#'
#' # Estimate the mean, 0.1 quantile, and 0.9 quantile with MCSEs using batch means.
#'
#' mcse(x)
#' mcse.q(x, 0.1)
#' mcse.q(x, 0.9)
#'
#' # Estimate the mean, 0.1 quantile, and 0.9 quantile with MCSEs using overlapping batch means.
#'
#' mcse(x, method = "obm")
#' mcse.q(x, 0.1, method = "obm")
#' mcse.q(x, 0.9, method = "obm")
#'
#' # Estimate E(x^2) with MCSE using spectral methods.
#'
#' g = function(x) { x^2 }
#' mcse(x, g = g, method = "tukey")
#'
#' @export

mcse = function(x, size = "sqroot", g = NULL, method = c("bm", "obm", "tukey", "bartlett"), warn = FALSE)
{
    if (! is.function(g))
        g = function(x) return(x)
    n = length(x)
    if (n < 1000)
    {
        if (warn)
            warning("too few samples (less than 1,000)")
        if (n < 10)
            return(NA)
    }
    if (size == "sqroot") 
    {
        b = floor(sqrt(n))
        a = floor(n / b)
    }
    else if (size == "cuberoot") 
    {
        b = floor(n^(1 / 3))
        a = floor(n / b)
    }
    else
    {
        if (! is.numeric(size) || size <= 1 || size == Inf)
            stop("'size' must be a finite numeric quantity larger than 1.")
        b = floor(size)
        a = floor(n / b)
    }
    method = match.arg(method)
    if (method == "bm")
    {
        y = sapply(1:a, function(k) return(mean(g(x[((k - 1) * b + 1):(k * b)]))))
        mu.hat = mean(g(x)) 
        var.hat = b * sum((y - mu.hat)^2) / (a - 1)
        se = sqrt(var.hat / n)
    }
    else if (method == "obm")
    {
        a = n - b + 1
        y = sapply(1:a, function(k) return(mean(g(x[k:(k + b - 1)]))))
        mu.hat = mean(g(x))
        var.hat = n * b * sum((y - mu.hat)^2) / (a - 1) / a
        se = sqrt(var.hat / n)
    } 
    else if (method == "tukey")
    {
        alpha = 1:b
        alpha = (1 + cos(pi * alpha / b)) / 2 * (1 - alpha / n)
        mu.hat = mean(g(x))
        R = sapply(0:b, function(j) return(mean((g(x[1:(n - j)]) - mu.hat) * (g(x[(j + 1):n]) - mu.hat))))
        var.hat = R[1] + 2 * sum(alpha * R[-1])
        se = sqrt(var.hat / n)
    }
    else # method == "bartlett"
    {
        alpha = 1:b
        alpha = (1 - abs(alpha) / b) * (1 - alpha / n)
        mu.hat = mean(g(x))
        R = sapply(0:b, function(j) return(mean((g(x[1:(n - j)]) - mu.hat) * (g(x[(j + 1):n]) - mu.hat))))
        var.hat = R[1] + 2 * sum(alpha * R[-1])
        se = sqrt(var.hat / n)
    }
    list(est = mu.hat, se = se)
}

#' Apply \code{mcse} to each column of a matrix or data frame of MCMC samples.
#'
#' @param x a matrix or data frame with each row being a draw from the multivariate distribution of interest.
#' @param size the batch size. The default value is \dQuote{\code{sqroot}}, which uses the square root of the sample size. \dQuote{\code{cuberoot}} will cause the function to use the cube root of the sample size. A numeric value may be provided if neither \dQuote{\code{sqroot}} nor \dQuote{\code{cuberoot}} is satisfactory.
#' @param g a function such that \eqn{E(g(x))} is the quantity of interest. The default is \code{NULL}, which causes the identity function to be used.
#' @param method the method used to compute the standard error. This is one of \dQuote{\code{bm}} (batch means, the default), \dQuote{\code{obm}} (overlapping batch means), \dQuote{\code{tukey}} (spectral variance method with a Tukey-Hanning window), or \dQuote{\code{bartlett}} (spectral variance method with a Bartlett window).
#' @return \code{mcse.mat} returns a matrix with \code{ncol(x)} rows and two columns. The row names of the matrix are the same as the column names of \code{x}. The column names of the matrix are \dQuote{\code{est}} and \dQuote{\code{se}}. The \eqn{j}th row of the matrix contains the result of applying \code{mcse} to the \eqn{j}th column of \code{x}.
#' @seealso
#' \code{\link{mcse}}, which acts on a vector.
#'
#' \code{\link{mcse.q}} and \code{\link{mcse.q.mat}}, which compute standard errors for quantiles.
#' @export

mcse.mat = function(x, size = "sqroot", g = NULL, method = c("bm", "obm", "tukey", "bartlett"))
{
    if (! is.matrix(x) && ! is.data.frame(x))
        stop("'x' must be a matrix or data frame.")
    num = ncol(x)
    vals = matrix(NA, num, 2)
    colnames(vals) = c("est", "se")
    rownames(vals) = colnames(x)
    res = apply(x, 2, mcse, size = size, g = g, method = method)
    for (i in 1:num)
        vals[i, ] = c(res[[i]]$est, res[[i]]$se)
    vals
}

#' Compute Monte Carlo standard errors for quantiles.
#'
#' @param x a vector of values from a Markov chain.
#' @param q the quantile of interest.
#' @param size the batch size. The default value is \dQuote{\code{sqroot}}, which uses the square root of the sample size. A numeric value may be provided if \dQuote{\code{sqroot}} is not satisfactory.
#' @param g a function such that the \eqn{q}th quantile of the univariate distribution function of \eqn{g(x)} is the quantity of interest. The default is \code{NULL}, which causes the identity function to be used.
#' @param method the method used to compute the standard error. This is one of \dQuote{\code{bm}} (batch means, the default), \dQuote{\code{obm}} (overlapping batch means), or \dQuote{\code{sub}} (subsampling bootstrap).
#' @param warn a logical value indicating whether the function should issue a warning if the sample size is too small (less than 1,000).
#' @return \code{mcse.q} returns a list with two elements:
#'         \item{est}{an estimate of the \eqn{q}th quantile of the univariate distribution function of \eqn{g(x)}.}
#'         \item{se}{the Monte Carlo standard error.}
#' @references
#' Flegal, J. M. (2012) Applicability of subsampling bootstrap methods in Markov chain Monte Carlo. In Wozniakowski, H. and Plaskota, L., editors, \emph{Monte Carlo and Quasi-Monte Carlo Methods 2010} (to appear). Springer-Verlag.
#'
#' Flegal, J. M. and Jones, G. L. (2010) Batch means and spectral variance estimators in Markov chain Monte Carlo. \emph{The Annals of Statistics}, \bold{38}, 1034--1070.
#'
#' Flegal, J. M. and Jones, G. L. (2011) Implementing Markov chain Monte Carlo: Estimating with confidence. In Brooks, S., Gelman, A., Jones, G. L., and Meng, X., editors, \emph{Handbook of Markov Chain Monte Carlo}, pages 175--197. Chapman & Hall/CRC Press.
#'
#' Flegal, J. M., Jones, G. L., and Neath, R. (2012) Markov chain Monte Carlo estimation of quantiles. \emph{University of California, Riverside, Technical Report}.
#'
#' Jones, G. L., Haran, M., Caffo, B. S. and Neath, R. (2006) Fixed-width output analysis for Markov chain Monte Carlo. \emph{Journal of the American Statistical Association}, \bold{101}, 1537--1547.
#' @seealso
#' \code{\link{mcse.q.mat}}, which applies \code{mcse.q} to each column of a matrix or data frame.
#'
#' \code{\link{mcse}} and \code{\link{mcse.mat}}, which compute standard errors for expectations.
#' @examples
#'
#' # Create 10,000 iterations of an AR(1) Markov chain with rho = 0.9.
#'
#' n = 10000
#' x = double(n)
#' x[1] = 2
#' for (i in 1:(n - 1))
#'     x[i + 1] = 0.9 * x[i] + rnorm(1)
#'
#' # Estimate the mean, 0.1 quantile, and 0.9 quantile with MCSEs using batch means.
#'
#' mcse(x)
#' mcse.q(x, 0.1)
#' mcse.q(x, 0.9)
#'
#' # Estimate the mean, 0.1 quantile, and 0.9 quantile with MCSEs using overlapping batch means.
#'
#' mcse(x, method = "obm")
#' mcse.q(x, 0.1, method = "obm")
#' mcse.q(x, 0.9, method = "obm")
#'
#' # Estimate E(x^2) with MCSE using spectral methods.
#'
#' g = function(x) { x^2 }
#' mcse(x, g = g, method = "tukey")
#'
#' @export

mcse.q = function(x, q, size = "sqroot", g = NULL, method = c("bm", "obm", "sub"), warn = FALSE)
{
    if (! is.function(g))
        g = function(x) return(x)
    n = length(x)
    if (n < 1000)
    {
        if (warn)
            warning("too few samples (less than 1,000)")
        if (n < 10)
            return(NA)
    }
    if (size == "sqroot") 
    {
        b = floor(sqrt(n))
        a = floor(n / b)
    }
    else
    {
        if (! is.numeric(size) || size <= 1 || size == Inf)
            stop("'size' must be a finite numeric quantity larger than 1.")
        b = floor(size)
        a = floor(n / b)
    }
    method = match.arg(method)
    counting = function(var.vector, var.number)
    {
        return(length(var.vector[var.vector <= var.number]))
    }
    if (! is.numeric(q) || q <= 0 || q >= 1)
        stop("'q' must be from (0, 1).")
    quant = function(input) { quantile(input, prob = q, type = 1, names = FALSE) }
    if (method == "bm")
    {
        xi.hat = quant(g(x))
        y = sapply(1:a, function(k) return(counting(g(x[((k - 1) * b + 1):(k * b)]), xi.hat))) / b
        mu.hat = mean(y)
        var.hat = b * sum((y - mu.hat)^2) / (a - 1)
        f.hat.junk = density(g(x), from = xi.hat, to = xi.hat, n = 1)
        f.hat = f.hat.junk$y
        se = sqrt(var.hat / n) / f.hat
    }
    else if (method == "obm")
    {
        xi.hat = quant(g(x))
        a = n - b + 1
        y = sapply(1:a, function(k) return(counting(g(x[k:(k + b - 1)]), xi.hat))) / b
        mu.hat = mean(y)
        var.hat = n * b * sum((y - mu.hat)^2) / (a - 1) / a
        f.hat.junk = density(g(x), from = xi.hat, to = xi.hat, n = 1)
        f.hat = f.hat.junk$y
        se = sqrt(var.hat / n) / f.hat
    } 
    else # method == "sub"
    {
        xi.hat = quant(g(x))
        a = n - b + 1
        y = sapply(1:a, function(k) return(quant(g(x[k:(k + b - 1)]))))
        mu.hat = mean(y)
        var.hat = n * b * sum((y - mu.hat)^2) / (a - 1) / a
        se = sqrt(var.hat / n)
    }
    list(est = xi.hat, se = se)      
}

#' Apply \code{mcse.q} to each column of a matrix or data frame of MCMC samples.
#'
#' @param x a matrix or data frame with each row being a draw from the multivariate distribution of interest.
#' @param q the quantile of interest.
#' @param size the batch size. The default value is \dQuote{\code{sqroot}}, which uses the square root of the sample size. \dQuote{\code{cuberoot}} will cause the function to use the cube root of the sample size. A numeric value may be provided if \dQuote{\code{sqroot}} is not satisfactory.
#' @param g a function such that the \eqn{q}th quantile of the univariate distribution function of \eqn{g(x)} is the quantity of interest. The default is \code{NULL}, which causes the identity function to be used.
#' @param method the method used to compute the standard error. This is one of \dQuote{\code{bm}} (batch means, the default), \dQuote{\code{obm}} (overlapping batch means), or \dQuote{\code{sub}} (subsampling bootstrap).
#' @return \code{mcse.q.mat} returns a matrix with \code{ncol(x)} rows and two columns. The row names of the matrix are the same as the column names of \code{x}. The column names of the matrix are \dQuote{\code{est}} and \dQuote{\code{se}}. The \eqn{j}th row of the matrix contains the result of applying \code{mcse.q} to the \eqn{j}th column of \code{x}.
#' @seealso \code{\link{mcse.q}}, which acts on a vector.
#'
#' \code{\link{mcse}} and \code{\link{mcse.mat}}, which compute standard errors for expectations.
#' @export

mcse.q.mat = function(x, q, size = "sqroot", g = NULL, method = c("bm", "obm", "sub"))
{
    if (! is.matrix(x) && ! is.data.frame(x))
        stop("'x' must be a matrix or data frame.")
    num = ncol(x)
    vals = matrix(NA, num, 2)
    colnames(vals) = c("est", "se")
    rownames(vals) = colnames(x)
    res = apply(x, 2, mcse.q, q = q, size = size, g = g, method = method)
    for (i in 1:num)
        vals[i, ] = c(res[[i]]$est, res[[i]]$se)
    vals
}

#' Create a plot that shows how Monte Carlo estimates change with increasing sample size.
#'
#' @param x a sample vector.
#' @param g a function such that \eqn{E(g(x))} is the quantity of interest. The default is \code{g = \link{mean}}.
#' @param main an overall title for the plot. The default is \dQuote{\code{Estimates vs Sample Size}}.
#' @param add logical. If \code{TRUE}, add to a current plot.
#' @param \dots additional arguments to the plotting function.
#' @return \code{NULL}
#' @examples
#' \dontrun{
#' estvssamp(x, main = expression(E(beta)))
#' estvssamp(y, add = TRUE, lty = 2, col = "red")}
#' @export

estvssamp = function(x, g = mean, main = "Estimates vs Sample Size", add = FALSE,...)
{
    if (length(x) < 100)
        size = 1
    else
        size = length(x) %/% 100
    n = seq(size, length(x), by = size)
    est = c()
    for (j in n)
        est = c(est, g(x[1:j]))
    if (add)
        lines(n, est,...)
    else
        plot(n, est, main = main, type = "l", xlab = "Sample Size", ylab = "MC Estimate",...)
}

#' Estimate effective sample size (ESS) as described in Kass et al. (1998) and Robert and Casella (2004; p. 500).
#'
#' @details ESS is the size of an iid sample with the same variance as the current sample. ESS is given by \deqn{\mbox{ESS}=T/\eta,}{ESS = T / \eta,} where \deqn{\eta=1+2\sum \mbox{all lag autocorrelations}.}{\eta = 1 + 2 \sum all lag autocorrelations.}
#'
#' @param x a vector of values from a Markov chain.
#' @param imse logical. If \code{TRUE}, use an approach that is analogous to Geyer's initial monotone positive sequence estimator (IMSE), where correlations beyond a certain lag are removed to reduce noise.
#' @param verbose logical. If \code{TRUE} and \code{imse = TRUE}, inform about the lag at which truncation occurs, and warn if the lag is probably too small. 
#' @return The function returns the estimated effective sample size.
#' @references
#' Kass, R. E., Carlin, B. P., Gelman, A., and Neal, R. (1998) Markov chain Monte Carlo in practice: A roundtable discussion. \emph{The American Statistician}, \bold{52}, 93--100.
#'
#' Robert, C. P. and Casella, G. (2004) \emph{Monte Carlo Statistical Methods}. New York: Springer.
#'
#' Geyer, C. J. (1992) Practical Markov chain Monte Carlo. \emph{Statistical Science}, \bold{7}, 473--483.
#' @export

# ess = function(x, imse = TRUE, verbose = FALSE)
# {
#     if (imse)
#     {
#         chain.acov = acf(x, type = "covariance", plot = FALSE)$acf
#         len = length(chain.acov)
#         gamma.acov = chain.acov[1:(len - 1)] + chain.acov[2:len]
#         k = 1
#         while (k < length(gamma.acov) && gamma.acov[k + 1] > 0 && gamma.acov[k] >= gamma.acov[k + 1])
#             k = k + 1
#         if (verbose)
#             cat("truncated after ", k, " lags\n")
#         if (k == length(gamma.acov) && verbose)
#             warning("may need to compute more autocovariances/autocorrelations for ess")
#         if (k == 1)
#             time = 1
#         else
#         {
#             chain.acor = acf(x, type = "correlation", plot = FALSE)$acf
#             time = 1 + 2 * sum(chain.acor[2:k])
#         }
#     }
#     else
#     {
#         chain.acor = acf(x, type = "correlation", plot = FALSE)$acf
#         time = 1 + 2 * sum(chain.acor[-1])
#     }
#     length(x) / time
# }

ess <- function(x, g = NULL)
{
    chain <- as.matrix(x)

    if (is.function(g)) 
        chain <- t(apply(x, 1, g))

    n <- dim(chain)[1]
    p <- dim(chain)[2]

    lambda <- apply(chain, 2, var)

    sigma <- (mcse.mat(chain)[,2])^2*n

    return(n*lambda/sigma)
}

