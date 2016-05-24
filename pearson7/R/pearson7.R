
#' Evaluate the density for the Pearson VII distribution with shape parameter 3/2.
#'
#' @details If \code{mu} is not specified, it assumes the default value of 0. If \code{sigma} is not specified, it assumes the default value of 1.
#'
#' The Pearson VII distribution with location \eqn{\mu}, scale \eqn{\sigma}, and shape 3/2 has density \deqn{f(x)=1/(2\sigma)[1+\{(x-\mu)/\sigma\}^2]^{-3/2}.}{f(x)=1/(2\sigma)[1+{(x-\mu)/\sigma}^2]^{-3/2}.}
#' @param x vector of quantiles.
#' @param mu vector of means.
#' @param sigma vector of scales.
#' @param log logical; if \code{TRUE}, probabilities p are given as log(p).
#' @return the density.
#' @references
#' Hughes, J., Shastry, S., Hancock, W. O., and Fricks, J. (2013) Estimating velocity for processive motor proteins with random detachment. \emph{Journal of Agricultural, Biological, and Environmental Statistics}, in press.
#' @references
#' Pearson, K. (1916) Mathematical contributions to the theory of evolution. xix. second supplement to a memoir on skew variation. \emph{Philosophical Transactions of the Royal Society of London. Series A, Containing Papers of a Mathematical or Physical Character}, \bold{216}, 429--457.
#' @seealso \code{\link{ppearson7}}, \code{\link{qpearson7}}, \code{\link{rpearson7}}
#' @examples
#'
#' curve(dpearson7(x), -5, 5, lwd = 2, n = 500, ylab = "f(x)")
#' curve(dnorm(x), lwd = 2, lty = 2, n = 500, add = TRUE)
#'
#' @export

dpearson7 = function(x, mu = 0, sigma = 1, log = FALSE)
{
    if (log)
        f = -1.5 * log(1 + ((x - mu) / sigma)^2) - log(sigma) - log(2)
    else
        f = 1 / (sigma * 2) * (1 + ((x - mu) / sigma)^2)^(-1.5)
    f
}

#' Evaluate the distribution function for the Pearson VII distribution with shape parameter 3/2.
#'
#' @details If \code{mu} is not specified, it assumes the default value of 0. If \code{sigma} is not specified, it assumes the default value of 1.
#'
#' The Pearson VII distribution with location \eqn{\mu}, scale \eqn{\sigma}, and shape 3/2 has cdf \deqn{F(x)=\{1+(x-\mu)/\sqrt{\sigma^2+(x-\mu)^2}\}/2.}{F(x)=[1+(x-\mu)/\sqrt{\sigma^2+(x-\mu)^2}]/2.}
#' @param q vector of quantiles.
#' @param mu vector of means.
#' @param sigma vector of scales.
#' @param log.p logical; if \code{TRUE}, probabilities p are given as log(p).
#' @param lower.tail logical; if \code{TRUE} (default), probabilities are \eqn{P(X\le x)}, otherwise \eqn{P(X>x)}.
#' @return the probability.
#' @references
#' Hughes, J., Shastry, S., Hancock, W. O., and Fricks, J. (2013) Estimating velocity for processive motor proteins with random detachment. \emph{Journal of Agricultural, Biological, and Environmental Statistics}, in press.
#' @seealso \code{\link{dpearson7}}, \code{\link{qpearson7}}, \code{\link{rpearson7}}
#' @examples
#'
#' curve(ppearson7(x), 0, 5, lwd = 2, ylim = c(0.8, 1), ylab = "F(x)")
#' curve(pnorm(x), lwd = 2, lty = 2, add = TRUE)
#'
#' @export

ppearson7 = function(q, mu = 0, sigma = 1, lower.tail = TRUE, log.p = FALSE)
{
    f = (1 + (q - mu) / (sqrt(sigma^2 + (q - mu)^2))) / 2
    if (! lower.tail)
        f = 1 - f
    if (log.p)
        f = log(f)
    f
}

#' Evaluate the quantile function for the Pearson VII distribution with shape parameter 3/2.
#'
#' @details If \code{mu} is not specified, it assumes the default value of 0. If \code{sigma} is not specified, it assumes the default value of 1.
#'
#' The Pearson VII distribution with location \eqn{\mu}, scale \eqn{\sigma}, and shape 3/2 has quantile function \deqn{F^{-1}(x)=\mu+(\sigma/2)(2x-1)/\sqrt{x(1-x)}.}
#' @param p vector of probabilities.
#' @param mu vector of means.
#' @param sigma vector of scales.
#' @param log.p logical; if \code{TRUE}, probabilities p are given as log(p).
#' @param lower.tail logical; if \code{TRUE} (default), probabilities are \eqn{P(X\le x)}, otherwise \eqn{P(X>x)}.
#' @return the quantile.
#' @references
#' Hughes, J., Shastry, S., Hancock, W. O., and Fricks, J. (2013) Estimating velocity for processive motor proteins with random detachment. \emph{Journal of Agricultural, Biological, and Environmental Statistics}, in press.
#' @seealso \code{\link{dpearson7}}, \code{\link{ppearson7}}, \code{\link{rpearson7}}
#' @examples
#'
#' curve(qpearson7(x), 0, 1, lwd = 2, ylab = expression(F^{-1}*(x)))
#' curve(qnorm(x), lwd = 2, lty = 2, n = 500, add = TRUE)
#'
#' @export

qpearson7 = function(p, mu = 0, sigma = 1, lower.tail = TRUE, log.p = FALSE)
{
	if (! lower.tail)
	    p = 1 - p
	if (log.p)
	    p = exp(p)
    z = 2 * p - 1
    f = mu + sigma * z / sqrt(1 - z^2)
    f
}

#' Generate random deviates from a Pearson VII distribution with shape parameter 3/2.
#'
#' @details If \code{mu} is not specified, it assumes the default value of 0. If \code{sigma} is not specified, it assumes the default value of 1.
#' @param n number of observations.
#' @param mu vector of means.
#' @param sigma vector of scales.
#' @return random deviates.
#' @references
#' Hughes, J., Shastry, S., Hancock, W. O., and Fricks, J. (2013) Estimating velocity for processive motor proteins with random detachment. \emph{Journal of Agricultural, Biological, and Environmental Statistics}, in press.
#' @references
#' Devroye, L. (1986) \emph{Non-Uniform Random Variate Generation}. New York: Springer-Verlag.
#' @seealso \code{\link{dpearson7}}, \code{\link{ppearson7}}, \code{\link{qpearson7}}
#' @examples
#'
#' y = rpearson7(1000)
#' hist(y, prob = TRUE, breaks = 100, col = "gray")
#' curve(dpearson7(x), lwd = 2, col = "blue", add = TRUE)
#'
#' @export

rpearson7 = function(n, mu = 0, sigma = 1)
{
    sigma * rnorm(n) / sqrt(2 * rexp(n)) + mu
}

#' Compute the negative log likelihood for a sample.
#'
#' @details This function computes the negative log likelihood for \eqn{(\mu,\sigma)} given a sample. This function can be optimized using \code{\link{optim}}, but it is better to use \code{\link{pearson7.fit}}.
#' @param params a vector of parameter values.
#' @param y a vector of observations.
#' @return the negative log likelihood.
#' @seealso \code{\link{dpearson7}}, \code{\link{pearson7.fit}}
#' @export

pearson7.objective = function(params, y)
{
    mu = params[1]
    sigma = params[2]
    -sum(dpearson7(y, mu, sigma, log = TRUE))
}

# This function computes the Hessian.

pearson7.hessian = function(y, mu, sigma)
{
    n = length(y)
    z = y - mu
    h11 = 3 * sum((z^2 - sigma^2) / (sigma^2 + z^2)^2)
    h12 = -6 * sigma * sum(z / (sigma^2 + z^2)^2)
    h22 = n / sigma^2 - 9 * sum(z^2 / (sigma^2 + z^2)^2) - 3 / sigma^2 * sum(z^4 / (sigma^2 + z^2)^2)
    matrix(c(h11, h12, h12, h22), 2, 2, byrow = TRUE)
}

# This function completes one step of the Newton-Raphson algorithm.

pearson7.step = function(y, mu0, sigma0)
{
    delta = c(mu0, sigma0)
    n = length(y)
    z = y - mu0
    grad = c(3 * sum(z / (sigma0^2 + z^2)), -n / sigma0 + 3 / sigma0 * sum(z^2 / (sigma0^2 + z^2)))
    H = pearson7.hessian(y, mu0, sigma0)
    ddelta = solve(H) %*% grad
    delta - ddelta
}

#' Find the MLE for a sample from the Pearson VII distribution with shape parameter 3/2.
#'
#' @details This function uses a Newton-Raphson algorithm to find the MLE. The starting values for \eqn{\mu} and \eqn{\sigma} are the sample median and \eqn{\sqrt{3}}{\sqrt3} times the sample MAD, respectively. See the reference for details.
#' @param y a vector of observations.
#' @param mu0 an initial value for \eqn{\mu}.
#' @param sigma0 an initial value for \eqn{\sigma}.
#' @param tol the convergence tolerance.
#' @return \code{pearson7.fit} returns an object of class \dQuote{\code{pearson7}}, which is a list containing the following components.
#'         \item{theta.hat}{the estimates of \eqn{\mu} and \eqn{\sigma}.}
#'         \item{hessian}{the Hessian matrix evaluated at \code{theta.hat}.}
#'         \item{iterations}{the number of iterations required to attain convergence.}
#'         \item{value}{the value of the log likelihood at \code{theta.hat}.}
#' @references
#' Hughes, J., Shastry, S., Hancock, W. O., and Fricks, J. (2013) Estimating velocity for processive motor proteins with random detachment. \emph{Journal of Agricultural, Biological, and Environmental Statistics}, in press.
#' @seealso \code{\link{pearson7.objective}}
#' @examples
#'
#' y = rpearson7(100, 100, 10)
#' fit = pearson7.fit(y)
#' fit
#' summary(fit)
#'
#' @export

pearson7.fit = function(y, mu0 = median(y), sigma0 = sqrt(3) * median(abs(y - median(y))), tol = 1e-8)
{
    count = 2
    current = pearson7.step(y, mu0, sigma0)
    nxt = pearson7.step(y, current[1], current[2])
    while (max(abs(nxt - current)) > tol)
    {
        current = nxt
        nxt = pearson7.step(y, current[1], current[2])
        count = count + 1
    }
    nxt = as.vector(nxt)
    names(nxt) = c("mu", "sigma")
    fit = list(theta.hat = nxt, hessian = -pearson7.hessian(y, nxt[1], nxt[2]),
               iterations = count, value = -pearson7.objective(nxt, y))
    class(fit) = "pearson7"
    fit
}

#' Print a summary of a Pearson VII fit.
#'
#' @details This function displays (1) a table of estimates, (2) the value of the log likelihood, and (3) the number of Newton-Raphson iterations. Each row of the table of estimates shows the parameter estimate and the approximate \eqn{(1-\alpha)100\%}{(1-\alpha)100\%} confidence interval for the parameter.
#' @param object an object of class \dQuote{\code{pearson7}}, the result of a call to \code{\link{pearson7.fit}}.
#' @param alpha the significance level used to compute the confidence intervals. The default is 0.05.
#' @param digits the number of significant digits to display. The default is 4.
#' @param \dots additional arguments.
#' @seealso \code{\link{pearson7.fit}}
#' @method summary pearson7
#' @export

summary.pearson7 = function(object, alpha = 0.05, digits = 4, ...)
{
    ci = matrix(, 2, 2)
    se = sqrt(diag(solve(object$hessian)))
    z = qnorm(1 - alpha / 2)
    ci[1, ] = c(object$theta.hat[1] - z * se[1], object$theta.hat[1] + z * se[1])
    ci[2, ] = c(object$theta.hat[2] - z * se[2], object$theta.hat[2] + z * se[2])
    coef.table = cbind(object$theta.hat, ci)
    colnames(coef.table) = c("Estimate", "Lower", "Upper")
    rownames(coef.table) = names(object$theta.hat)
    cat("\nEstimates:\n")
    print(signif(coef.table, digits))
    cat("\nLog likelihood:", signif(object$value, digits), "\n\nNumber of iterations:", object$iterations, "\n\n")
}


