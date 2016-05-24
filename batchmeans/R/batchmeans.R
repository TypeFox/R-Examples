
#' Perform consistent batch means estimation on a vector of values from a Markov chain.
#'
#' @param x a vector of values from a Markov chain.
#' @param size the batch size. The default value is \dQuote{\code{sqroot}}, which uses the square root of the sample size. \dQuote{\code{cuberoot}} will cause the function to use the cube root of the sample size. A numeric value may be provided if neither \dQuote{\code{sqroot}} nor \dQuote{\code{cuberoot}} is satisfactory.
#' @param warn a logical value indicating whether the function should issue a warning if the sample size is too small (less than 1,000).
#' @return \code{bm} returns a list with two elements:
#'         \item{est}{the mean of the vector.}
#'         \item{se}{the MCMC standard error based on the consistent batch means estimator.}
#' @references
#' Jones, G. L., Haran, M., Caffo, B. S. and Neath, R. (2006) Fixed-width output analysis for Markov chain Monte Carlo. \emph{Journal of the American Statistical Association}, \bold{101}, 1537--1547.
#'
#' The following article is less technical and contains a direct comparison to the Gelman-Rubin diagnostic.
#'
#' Flegal, J. M., Haran, M. and Jones, G. L. (2008) Markov chain Monte Carlo: Can we trust the third significant figure? \emph{Statistical Science}, \bold{23}, 250--260.
#' @seealso \code{\link{bmmat}}, which applies \code{bm} to each column of a matrix or data frame.
#' @export

bm = function(x, size = "sqroot", warn = FALSE)
{
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
    y = sapply(1:a, function(k) return(mean(x[((k - 1) * b + 1):(k * b)])))
    mu.hat = mean(y)
    var.hat = b * sum((y - mu.hat)^2) / (a - 1)
    se = sqrt(var.hat / n)
    list(est = mu.hat, se = se)
}

#' Apply \code{bm} to each column of a matrix or data frame of MCMC samples.
#'
#' @param x a matrix or data frame with each row being a draw from the multivariate distribution of interest.
#' @return \code{bmmat} returns a matrix with \code{ncol(x)} rows and two columns. The row names of the matrix are the same as the column names of \code{x}. The column names of the matrix are \dQuote{\code{est}} and \dQuote{\code{se}}. The \eqn{j}th row of the matrix contains the result of applying \code{bm} to the \eqn{j}th column of \code{x}.
#' @seealso \code{\link{bm}}, which performs consistent batch means estimation for a vector.
#' @export

bmmat = function(x)
{
    if (! is.matrix(x) && ! is.data.frame(x))
        stop("'x' must be a matrix or data frame.")
    num = ncol(x)
    bmvals = matrix(NA, num, 2)
    colnames(bmvals) = c("est", "se")
    rownames(bmvals) = colnames(x)
    bmres = apply(x, 2, bm)
    for (i in 1:num)
        bmvals[i, ] = c(bmres[[i]]$est, bmres[[i]]$se)
    bmvals
}

#' Create a plot that shows how Monte Carlo estimates change with increasing sample size.
#'
#' @param x a sample vector.
#' @param fun a function such that \eqn{E(fun(x))} is the quantity of interest. The default is \code{fun = \link{mean}}.
#' @param main an overall title for the plot. The default is \dQuote{\code{Estimates vs Sample Size}}.
#' @param add logical. If \code{TRUE}, add to a current plot.
#' @param \dots additional arguments to the plotting function.
#' @return \code{NULL}
#' @examples
#' \dontrun{
#' estvssamp(x, main = expression(E(beta)))
#' estvssamp(y, add = TRUE, lty = 2, col = "red")}
#' @export

estvssamp = function(x, fun = mean, main = "Estimates vs Sample Size", add = FALSE,...)
{
    if (length(x) < 100)
        size = 1
    else
        size = length(x) %/% 100
    n = seq(size, length(x), by = size)
    est = c()
    for (j in n)
        est = c(est, fun(x[1:j]))
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

ess = function(x, imse = TRUE, verbose = FALSE)
{
    if (imse)
    {
        chain.acov = acf(x, type = "covariance", plot = FALSE)$acf
        len = length(chain.acov)
        gamma.acov = chain.acov[1:(len - 1)] + chain.acov[2:len]
        k = 1
        while (k < length(gamma.acov) && gamma.acov[k + 1] > 0 && gamma.acov[k] >= gamma.acov[k + 1])
            k = k + 1
        if (verbose)
            cat("truncated after ", k, " lags\n")
        if (k == length(gamma.acov) && verbose)
            warning("may need to compute more autocovariances/autocorrelations for ess")
        if (k == 1)
            time = 1
        else
        {
            chain.acor = acf(x, type = "correlation", plot = FALSE)$acf
            time = 1 + 2 * sum(chain.acor[2:k])
        }
    }
    else
    {
        chain.acor = acf(x, type = "correlation", plot = FALSE)$acf
        time = 1 + 2 * sum(chain.acor[-1])
    }
    length(x) / time
}




