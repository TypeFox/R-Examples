
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA


################################################################################
# FUNCTION:             DESCRIPTION:
#  hngarchSim            Simulates an HN-GARCH(1,1) Time Series Process
#  hngarchFit            Fits a HN-GARCH model by Gaussian Maximum Likelihood
#  print.hngarch         Print method, reports results
#  summary.hngarch       Summary method, diagnostic analysis
#  hngarchStats          Computes Unconditional Moments of a HN-GARCH Process
################################################################################



hngarchSim =
function(model = list(lambda = 4, omega = 4*0.0002, alpha = 0.3*0.0002,
beta = 0.3, gamma = 0, rf = 0), n = 1000, innov = NULL, n.start = 100,
start.innov = NULL, rand.gen = rnorm, ...)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Simulates a HN-GARCH time series with user supplied innovations.

    # Details:
    #   The function simulates a Heston Nandi Garch(1,1) process with
    #   structure parameters specified through the list
    #   `model(lambda, omega, alpha, beta, gamma, rf)'
    #   The function returns the simulated time series points
    #   neglecting those from the first "start.innov" innovations.

    # Example:
    #   x = hngarch()
    #   plot(100*x, type="l", xlab="Day numbers",
    #     ylab="Daily Returns %", main="Heston Nandi GARCH")
    #   S0 = 1
    #   plot(S0*exp(cumsum(x)), type="l", xlab="Day Numbers",
    #     ylab="Daily Prices", main="Heston Nandi GARCH") }

    # FUNCTION:

    # Innovations:
    if (is.null(innov)) innov = rand.gen(n, ...)
    if (is.null(start.innov)) start.innov = rand.gen(n.start, ...)

    # Parameters:
    lambda = model$lambda
    omega = model$omega
    alpha = model$alpha
    beta = model$beta
    gamma = model$gamma
    rfr = model$rf

    # Start values:
    x = h = Z = c(start.innov, innov)
    nt = n.start + n

    # Recursion:
    h[1] = ( omega + alpha )/( 1 - alpha*gamma*gamma - beta )
    x[1] = rfr + lambda*h[1] + sqrt(h[1]) * Z[1]
    for (i in 2:nt) {
        h[i] = omega + alpha*(Z[i-1] - gamma*sqrt(h[i-1]))^2 + beta*h[i-1]
        x[i] = rfr + lambda*h[i] + sqrt(h[i]) * Z[i] }

    # Series:
    x = x[-(1:n.start)]

    # Return Value:
    x
}


# ------------------------------------------------------------------------------


hngarchFit =
function(x, model = list(lambda = -0.5, omega = var(x), alpha = 0.1*var(x),
beta = 0.1, gamma = 0, rf = 0), symmetric = TRUE, trace = FALSE, title =
NULL, description = NULL, ...)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Fits Heston-Nandi Garch(1,1) time series model

    # FUNCTION:

    # Parameters:
    rfr = model$rf
    lambda = model$lambda
    omega = model$omega
    alpha = model$alpha
    beta = model$beta
    gam = model$gamma

    # Continue:
    params = c(lambda = lambda, omega = omega, alpha = alpha,
        beta = beta, gamma = gam, rf = rfr)

    # Transform Parameters and Calculate Start Parameters:
    par.omega = -log((1-omega)/omega)  # for 2
    par.alpha = -log((1-alpha)/alpha)  # for 3
    par.beta = -log((1-beta)/beta)     # for 4
    par.start = c(lambda, par.omega, par.alpha, par.beta)
    if (!symmetric) par.start = c(par.start, gam)

    # Initial Log Likelihood:
    opt = list()
    opt$value = .llhHNGarch(par = par.start,
        trace = trace, symmetric = symmetric, rfr = rfr, x = x)
    opt$estimate = par.start
    if (trace) {
        print(c(lambda, omega, alpha, beta, gam))
        print(opt$value)
    }

    # Estimate Parameters:
    opt = nlm(.llhHNGarch, par.start,
        trace = trace, symmetric = symmetric, rfr = rfr, x = x, ...)

    # Log-Likelihood:
    opt$minimum = -opt$minimum + length(x)*sqrt(2*pi)
    opt$params = params
    opt$symmetric = symmetric

    # LLH, h, and z for Final Estimates:
    final = .llhHNGarch(opt$estimate, trace = FALSE, symmetric, rfr, x)
    opt$h = attr(final, "h")
    opt$z = attr(final, "Z")

    # Backtransform Estimated parameters:
    lambda = opt$estimate[1]
    omega = opt$estimate[2] = (1 / (1+exp(-opt$estimate[2])))
    alpha = opt$estimate[3] = (1 / (1+exp(-opt$estimate[3])))
    beta = opt$estimate[4] = (1 / (1+exp(-opt$estimate[4])))
    if (symmetric) opt$estimate[5] = 0
    gam = opt$estimate[5]
    names(opt$estimate) = c("lambda", "omega", "alpha", "beta", "gamma")

    # Add to Output:
    opt$model = list(lambda = lambda, omega = omega, alpha = alpha,
        beta = beta, gamma = gam, rf = rfr)
    opt$x = x

    # Statistics - Printing:
    opt$persistence = beta + alpha*gam*gam
    opt$sigma2 = ( omega + alpha ) / ( 1 - opt$persistence )

    # Print Estimated Parameters:
    if (trace) print(opt$estimate)

    # Call:
    opt$call = match.call()

    # Add title and description:
    if (is.null(title))
        title = "Heston-Nandi Garch Parameter Estimation"
    opt$title = title
    if (is.null(description))
        description = description()
    opt$description = description

    # Return Value:
    class(opt) = "hngarch"
    invisible(opt)
}


# ------------------------------------------------------------------------------


.llhHNGarch =
function(par, trace, symmetric, rfr, x)
{
    # h = sigma^2
    h = Z = x
    lambda = par[1]

    # Transform - to keep them between 0 and 1:
    omega = 1 / (1+exp(-par[2]))
    alpha = 1 / (1+exp(-par[3]))
    beta = 1 / (1+exp(-par[4]))

    # Add gamma if selected:
    if (!symmetric) gam = par[5] else gam = 0

    # HN Garch Filter:
    h[1] = ( omega + alpha )/( 1 - alpha*gam*gam - beta)
    Z[1] = ( x[1] - rfr - lambda*h[1] ) / sqrt(h[1])
    for ( i in 2:length(Z) ) {
        h[i] = omega + alpha * ( Z[i-1] - gam * sqrt(h[i-1]) )^2 +
            beta * h[i-1]
        Z[i] = ( x[i] - rfr - lambda*h[i] ) / sqrt(h[i])
    }

    # Calculate Log - Likelihood for Normal Distribution:
    llhHNGarch = -sum(log( dnorm(Z)/sqrt(h) ))
    if (trace) {
        cat("Parameter Estimate\n")
        print(c(lambda, omega, alpha, beta, gam))
    }

    # Attribute Z and h to the result:
    attr(llhHNGarch, "Z") = Z
    attr(llhHNGarch, "h") = h

    # Return Value:
    llhHNGarch
}


# ------------------------------------------------------------------------------


print.hngarch =
function(x, ...)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Print method for the  HN-GARCH time series model.

    # Arguments:
    #   x - an object of class "hngarch" as returned by the
    #   function "hngarchFit"

    # FUNCTION:

    # Print:
    object = x
    if (!inherits(object, "hngarch"))
        stop("method is only for garch objects")

    # Title:
    cat("\nTitle:\n ")
    cat(object$title, "\n")

    # Call:
    cat("\nCall:\n ", deparse(object$call), "\n", sep = "")

    # Parameters:
    cat("\nParameters:\n")
    print(format(object$params, digits = 4, ...), print.gap = 2,
        quote = FALSE)

    # Coefficients:
    cat("\nCoefficients: lambda, omega, alpha, beta, gamma\n")
    print(format(object$estimate, digits = 4, ...), print.gap = 2,
        quote = FALSE)

    # Likelihood:
    cat("\nLog-Likelihood:\n ")
    cat(object$minimum, "\n")

    # Persisitence and Variance:
    cat("\nPersistence and Variance:\n ")
    cat(object$persistence, "\n ")
    cat(object$sigma2, "\n")

    # Description:
    cat("\nDescription:\n ")
    cat(object$description, "\n\n")

    # Return Value:
    invisible()
}


# ------------------------------------------------------------------------------


summary.hngarch =
function(object, ...)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Summary method,
    #   Computes diagnostics for a HN-GARCH time series model.

    # Arguments:
    #   object - an object of class "hngarch" as returned by the
    #   function "hngarchFit"

    # FUNCTION:

    # Print:
    if (!inherits(object, "hngarch"))
        stop("method is only for garch objects")

    # Title:
    cat("\nTitle:\n")
    cat(object$title, "\n")

    # Call:
    cat("\nCall:\n", deparse(object$call), "\n", sep = "")

    # Parameters:
    cat("\nParameters:\n")
    print(format(object$params, digits = 4, ...), print.gap = 2,
        quote = FALSE)

    # Coefficients:
    cat("\nCoefficients: lambda, omega, alpha, beta, gamma\n")
    print(format(object$estimate, digits = 4, ...), print.gap = 2,
        quote = FALSE)

    # Likelihood:
    cat("\nLog-Likelihood:\n")
    cat(object$minimum, "\n")

    # Persisitence and Variance:
    cat("\nPersistence and Variance:\n")
    cat(object$persistence, "\n")
    cat(object$sigma2, "\n")

    # Create Graphs:
    plot(x = object$x, type = "l", xlab = "Days", ylab = "log-Returns",
        main = "Log-Returns", ...)
    plot(sqrt(object$h), type = "l", xlab = "Days", ylab = "sqrt(h)",
        main = "Conditional Standard Deviations", ...)
    # ... there are not resiudal yet implemented:
    # plot(object$residuals, type = "l", xlab = "Days", ylab = "Z",
    #     main = "Residuals", ...)

    # Return Value:
    invisible()
}


################################################################################


hngarchStats =
function(model)
{   # A function implemented by Diethelm Wuertz

    # Description:

    # Details:
    #   Calculates the first 4 moments of the unconditional log
    #   return distribution for a stationary HN GARCH(1,1) process
    #   with standard normally distributed innovations. The function
    #   returns a list with the theoretical values for the mean, the
    #   variance, the skewness and the kurtosis} of the (unconditional)
    #   log return distribution. We have also access to the persistence
    #   of the corresponding HN GARCH(1,1) process and to the values
    #   for E[sigma^2], E[sigma^4], E[sigma^6], and E[sigma^8], which are
    #   needed for the computation of the moments of the unconditional
    #   log return distribution. The only arguments are the risk free
    #   interest rate r and the structure parameters of the HN GARCH(1,1)
    #   process, which are specified in the model list model=list(alpha,
    #   beta, omega, gamma, lambda)}.

    # Reference:
    #   A function originally written by Reto Angliker
    #   License: GPL

    # Arguments:
    #   model - a moel specification for a Heston-Nandi Garch
    #       process.

    # FUNCTION:

    # Check:
    if (model$alpha < 0) {
        warning("Negative value for the parameter alpha")}
    if (model$beta < 0)
        {warning("Negative value for the parameter beta") }
    if (model$omega < 0)
        {warning("Negative value for the parameter omega")}

    # Short:
    lambda = model$lambda
    omega = model$omega
    alpha = model$alpha
    beta = model$beta
    gamma = model$gamma

    # Moments of the Normal Distribution
    expect2 = 1
    expect4 = 3
    expect6 = 15
    expect8 = 105

    # Symmetric Case:
    if(model$gamma == 0) {
        persistence = beta
        meansigma2 = (omega+alpha) /(1-beta)
        meansigma4 = (omega^2 + 2*omega*alpha + 2*omega*beta*meansigma2 +
            3*alpha^2 + 2*alpha*beta*meansigma2) / (1 - beta^2)
        meansigma6 = (omega^3 + 3*omega^2*alpha + 3*omega^2*beta*meansigma2 +
            9*omega*alpha^2 + 6*omega*alpha*beta*meansigma2 +
            3*omega*beta^2*meansigma4 + 15*alpha^3 +
            9*alpha^2*beta*meansigma2 + 3*alpha*beta^2*meansigma4) / (1-beta^3)
        meansigma8 =
            (omega^4 + expect8*alpha^4 + 12*omega^2*alpha*beta*meansigma2 +
            60*alpha^3*beta*meansigma2 + 18*alpha^2*beta^2*meansigma4 +
            4*alpha*beta^3*meansigma6 + 36*omega*alpha^2*beta*meansigma2 +
            12*omega*alpha*beta^2*meansigma4 + 4*omega^3*alpha +
            4*omega^3*beta*meansigma2 + 18*omega^2*alpha^2 +
            6*omega^2*beta^2*meansigma4 + 60*omega*alpha^3 +
            4*omega*beta^3*meansigma6)/ (1 - beta^4) }

    # Asymmetric Case:
    if(gamma != 0) {
        persistence = beta + alpha*gamma^2
        meansigma2 = (omega+alpha) / (1-beta-alpha*gamma^2)
        meansigma4 = (omega^2 + 2*omega*beta*meansigma2 +
            alpha^2*expect4 + 2*beta*meansigma2*alpha*expect2 +
            6*alpha^2*expect2*gamma^2*meansigma2 +
            2*omega*alpha*gamma^2*meansigma2 + 2*omega*alpha*expect2) /
            (1 - beta^2 - 2*beta*alpha*gamma^2 - alpha^2*gamma^4)
        meansigma6 =
            (3*omega*alpha^2*expect4 + 3*omega^2*alpha*gamma^2*meansigma2 +
            3*beta*meansigma2*alpha^2*expect4 +
            3*beta^2*meansigma4*alpha*expect2 +
            15*alpha^3*expect4*gamma^2*meansigma2 +
            15*alpha^3*expect2*gamma^4*meansigma4 +
            3*omega*alpha^2*gamma^4*meansigma4 +
            3*omega^2*beta*meansigma2 + 3*omega^2*alpha*expect2 +
            3*omega*beta^2*meansigma4 + omega^3 + alpha^3*expect6 +
            18*beta*meansigma4*alpha^2*expect2*gamma^2 +
            6*omega*beta*meansigma2*alpha*expect2 +
            6*omega*beta*meansigma4*alpha*gamma^2 +
            18*omega*alpha^2*expect2*gamma^2*meansigma2) /
            (1 - 3*beta^2*alpha*gamma^2 - 3*beta*alpha^2*gamma^4 -
            alpha^3*gamma^6 - beta^3)
        meansigma8 = (omega^4 + alpha^4*expect8 +
            6*omega^2*alpha^2*expect4 + 4*omega^3*beta*meansigma2 +
            4*omega^3*alpha*expect2 + 6*omega^2*beta^2*meansigma4 +
            4*omega*beta^3*meansigma6 + 4*omega*alpha^3*expect6 +
            12*omega^2*beta*meansigma2*alpha*expect2 +
            12*omega^2*beta*meansigma4*alpha*gamma^2 +
            36*omega^2*alpha^2*expect2*gamma^2*meansigma2 +
            4*omega^3*alpha*gamma^2*meansigma2 +
            6*omega^2*alpha^2*gamma^4*meansigma4 +
            6*beta^2*meansigma4*alpha^2*expect4 +
            4*beta^3*meansigma6*alpha*expect2 +
            4*beta*meansigma2*alpha^3*expect6 +
            28*alpha^4*expect6*gamma^2*meansigma2 +
            70*alpha^4*expect4*gamma^4*meansigma4 +
            28*alpha^4*expect2*gamma^6*meansigma6 +
            4*omega*alpha^3*gamma^6*meansigma6 +
            60*beta*meansigma4*alpha^3*expect4*gamma^2 +
            60* beta*meansigma6*alpha^3*expect2*gamma^4 +
            36*beta^2*meansigma6*alpha^2*expect2*gamma^2 +
            12*omega*beta*meansigma2*alpha^2*expect4 +
            12*omega*beta^2*meansigma4*alpha*expect2 +
            12*omega*beta^2*meansigma6*alpha*gamma^2 +
            12*omega*beta*meansigma6*alpha^2*gamma^4 +
            60*omega*alpha^3*expect4*gamma^2*meansigma2 +
            60*omega*alpha^3*expect2*gamma^4*meansigma4 +
            72*omega*beta*meansigma4*alpha^2*expect2*gamma^2) /
            (1 - beta^4 - alpha^4*gamma^8 - 4*beta^3*alpha*gamma^2 -
            6*beta^2*alpha^2*gamma^4 - 4*beta*alpha^3*gamma^6 ) }
    if (persistence > 1) { warning(paste(
        "The selected HN GARCH model is not stationary and",
        "the expressions for the moments are no more valid")) }

    # Leverage:
    leverage = -2*alpha*gamma*meansigma2

    # Unconditional Values:
    uc.mean = lambda*meansigma2
    uc.variance = lambda^2*(meansigma4 - meansigma2^2) + meansigma2
    uc.skewness = (3*lambda*meansigma4 - 3*lambda*meansigma2^2 +
        lambda^3*meansigma6 - 3*lambda^3*meansigma2*meansigma4 +
        2*lambda^3*meansigma2^3 ) / sqrt(uc.variance)^3
    uc.kurtosis = (meansigma4*3 + 6*lambda^2*meansigma6 -
        12*lambda^2*meansigma2*meansigma4 + 6*lambda^2*meansigma2^3 +
        lambda^4*meansigma8 - 4*lambda^4*meansigma2*meansigma6 +
        6*lambda^4*meansigma2^2*meansigma4 -
        3*lambda^4*meansigma2^4 ) / uc.variance^2

    # Return Value:
    list(mean = uc.mean, variance = uc.variance, skewness = uc.skewness,
        kurtosis = uc.kurtosis, persistence = persistence, leverage = leverage,
        meansigma2 = meansigma2, meansigma4 = meansigma4, meansigma6 =
        meansigma6, meansigma8 = meansigma8)
}


################################################################################

