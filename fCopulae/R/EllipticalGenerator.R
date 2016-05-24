
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
# FUNCTION:                  UTILITY FUNCTIONS:
#  ellipticalList             Returns list of implemented Elliptical copulae
#  ellipticalParam            Sets default parameters for an elliptical copula
#  ellipticalRange            Returns the range of valid rho values
#  ellipticalCheck            Checks if rho is in the valid range
# FUNCTION:                  ELLIPTICAL GENERATOR AND RELATED FUNCTIONS:
#  gfunc                      Generator function for elliptical distributions
#  gfuncSlider                Slider for generator, density and probability
#  .pelliptical               Univariate elliptical distribution probability
#  .delliptical               Univariate elliptical distribution density
#  .qelliptical               Univariate elliptical distribution quantiles
#  .qlogistic                 Fast tabulated logistic quantile function
#  .qlogisticData             Table generator for logistic quantiles
#  .qlogisticTable            Table for logistic quantiles
################################################################################


################################################################################
# UTILITY FUNCTIONS:
#  ellipticalParam            Sets Default parameters for an elliptical copula
#  ellipticalList             Returns list of implemented Elliptical copulae
#  ellipticalRange            Returns the range of valid rho values
#  ellipticalCheck            Checks if rho is in the valid range


ellipticalList <- 
    function()
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Returns list of implemented elliptical copulae

    # Arguments:

    # FUNCTION:

    # Compose List:
    ans = c("norm", "cauchy", "t", "logistic", "laplace", "kotz", "epower")

    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


ellipticalParam <-
    function(type = ellipticalList())
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Sets default parameters for elliptical copulae

    # Arguments:
    #   type -  a character string denoting the type of distribution.
    #       This may be either
    #       "norm" for the normal distribution, or
    #       "cauchy" for the Cauchy distribution, or
    #       "t" for the Student-t distribution, or
    #       "logistic" for the logistic distribution, or
    #       "laplace" for the distribution, or
    #       "kotz" for the original Kotz distribution, or
    #       "epower" for the exponential power distribution

    # Value:
    #   returns a list with two elements, 'param' sets the parameters
    #       which may be a vector, 'range' the range with minimum and
    #       maximum values for each of the parameters.

    # Example:
    #   ellipticalParam("norm"); ellipticalParam("t")

    # FUNCTION:

    # Settings:
    type = match.arg(type)

    # Parameter Values:
    #         ("norm", "cauchy",  "t", "logistic", "laplace", "kotz", "epower")
    lower  = c( -1,      -1,   -1,         -1,         -1,     -1,       -1)
    upper  = c( +1,      +1,   +1,         +1,         +1,     +1,       +1)
    rho    = c(3/4,     3/4,  3/4,        3/4,        3/4,    3/4,      3/4)
    param1 = c( NA,      NA,    4,         NA,         NA,      1,        1)
    param2 = c( NA,      NA,   NA,         NA,         NA,     NA,        1)

    # Create Parameter List:
    ans = list(type = type)
    if (type == "norm") {
        ans$param = c(rho = rho[1])
        ans$range = c(lower = lower[1], upper = upper[1])
    }
    if (type == "cauchy") {
        ans$param = c(rho = rho[2])
        ans$range = c(lower = lower[2], upper = upper[2])
    }
    if (type == "t") {
        ans$param = c(rho = rho[3], nu = param1[3])
        ans$range = c(lower = lower[3], upper = upper[3])
    }
    if (type == "logistic") {
        ans$param = c(rho = rho[4])
        ans$range = c(lower = lower[4], upper = upper[4])
    }
    if (type == "laplace") {
        ans$param = c(rho = rho[5])
        ans$range = c(lower = lower[5], upper = upper[5])
    }
    if (type == "kotz") {
        ans$param = c(rho = rho[6], r = param1[6])
        ans$range = c(lower = lower[6], upper = upper[6])
    }
    if (type == "epower") {
        ans$param = c(rho = rho[7], r = param1[7], s = param2[7])
        ans$range = c(lower = lower[7], upper = upper[7])
    }

    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


ellipticalRange <- 
    function(type = ellipticalList())
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Returns the range of valid alpha values

    # Arguments:
    #   type -  a character string denoting the type of distribution.
    #       This may be either
    #       "norm" for the normal distribution, or
    #       "cauchy" for the Cauchy distribution, or
    #       "t" for the Student-t distribution, or
    #       "logistic" for the logistic distribution, or
    #       "laplace" for the distribution, or
    #       "kotz" for the original Kotz distribution, or
    #       "epower" for the exponential power distribution

    # Example:
    #   ellipticalRange("norm"); ellipticalRange("t")

    # FUNCTION:

    # Type:
    type = match.arg(type)

    # Range:
    ans = ellipticalParam(type)$range
    attr(ans, "control") <- type

    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


ellipticalCheck <- 
    function(rho = 0.75, param = NULL, type = ellipticalList())
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Checks if alpha is in the valid range

    # Arguments:
    #   rho - correlation coefficient
    #   param - currently not used
    #   type -  a character string denoting the type of distribution.
    #       This may be either
    #       "norm" for the normal distribution, or
    #       "cauchy" for the Cauchy distribution, or
    #       "t" for the Student-t distribution, or
    #       "logistic" for the logistic distribution, or
    #       "laplace" for the distribution, or
    #       "kotz" for the original Kotz distribution, or
    #       "epower" for the exponential power distribution

    # Example:
    #   ellipticalCheck(0.5, NULL, "norm")
    #   ellipticalCheck(1.5, NULL, "t")

    # FUNCTION:

    # Type:
    type = match.arg(type)

    # Range:
    range = as.vector(ellipticalRange(type))
    if (rho < range[1] | rho > range[2]) {
        print(c(rho = rho))
        print(c(range = range))
        stop("rho is out of range")
    }

    # Return Value:
    invisible()
}


################################################################################
# FUNCTION:                  ELLIPTICAL GENERATOR AND RELATED FUNCTIONS:
#  gfunc                      Generator function for elliptical distributions
#  gfuncSlider                Slider for generator, density and probability
#  .pelliptical               Univariate elliptical distribution probability
#  .delliptical               Univariate elliptical distribution density
#  .qelliptical               Univariate elliptical distribution quantiles
#  .qlogistic                 Fast tabulated logistic quantile function
#  .qlogisticData             Table generator for logistic quantiles
#  .qlogisticTable            Table for logistic quantiles


gfunc <- 
    function(x, param = NULL, type = ellipticalList())
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Generator function for elliptical distributions

    # Arguments:
    #   x -  a numeric vector
    #   param - NULL, a numeric value, or a numeric vector adding.
    #       additional parameters to the generator function.
    #   type -  a character string denoting the type of distribution.
    #       This may be either
    #       "norm" for the normal distribution, or
    #       "cauchy" for the Cauchy distribution, or
    #       "t" for the Student-t distribution, or
    #       "logistic" for the logistic distribution, or
    #       "laplace" for the distribution, or
    #       "kotz" for the original Kotz distribution, or
    #       "epower" for the exponential power distribution

    # Value:
    #   Returns a numeric vector "g(x)" for the generator computed at
    #   the x values taken from the input vector. If x is missing,
    #   the normalizing constant "lambda" will be returned.

    # FUNCTION:

    # Match Arguments:
    type = match.arg(type)

    # Handle Missing x:
    if (missing(x)) {
        x = NA
        output = "lambda"
    } else {
        output = "g"
    }

    # Get Type:
    type = type[1]

    # Get Parameters:
    # if (is.null(param)) param = ellipticalParam$param

    # Create Generator:
    if (type == "norm") {
        g = exp(-x/2)
        lambda = 1 / (2*pi)
        param = NULL
    }
    if (type == "cauchy") {
        g = ( 1 + x )^ (-3/2 )
        lambda = 1 / (2*pi)
        param = NULL
    }
    if (type == "t") {
        if (is.null(param)) {
            nu = 4
        } else {
            nu = param[[1]]
        }
        g = ( 1 + x/nu )^ ( -(nu+2)/2 )
        lambda = 1/(2*pi)
        param = c(nu = nu)
    }
    if (type == "logistic"){
        g = exp(-x/2)/(1+exp(-x/2))^2
        # lambda:
        # integrate(function(x) { exp(-x)/(1+exp(-x))^2}, 0, Inf,
        #   subdivision = 10000, rel.tol = .Machine$double.eps^0.8)
        # 0.5 with absolute error < 2.0e-13
        lambda = 1 / pi
        param = NULL
    }
    if (type == "laplace") { # or "double exponential"
        # epower - with r = 1, s = 1
        # g = exp(-r*(x/2)^s)
        # lambda = s * r^(1/s) / ( 2 * pi * gamma(1/s) )
        g = exp(-sqrt(x))
        lambda = 1/(2*pi)
        param = NULL
    }
    if (type == "kotz") {
        # epower - with s = 1
        if (is.null(param)) {
            r = 1
        } else {
            r = param
        }
        g = exp(-r*(x/2))
        lambda = r/(2*pi)
        param = c(r = r)
    }
    if (type == "epower") {
        if (is.null(param)) {
            r = 1
            s = 1
        } else {
            r = param[[1]]
            s = param[[2]]
        }
        g = exp(-r*(x/2)^s)
        lambda = s * r^(1/s) / ( 2 * pi * gamma(1/s) )
        param = c(r = r, s = s)
    }

    # Output:
    output = output[1]
    if (output == "g") {
        ans = g
    } else if (output == "lambda") {
        ans = lambda
    }

    # Add Control:
    if (output == "g") {
        attr(ans, "control") <-
            c(copula = "elliptical", type = type, lambda = as.character(lambda))
    } else if (output == "lambda") {
        if (is.null(param)) {
            attr(ans, "control") <-
                unlist(list(copula = "elliptical", type = type))
        } else {
            attr(ans, "control") <-
                unlist(list(copula = "elliptical", type = type, param = param))
        }
    }

    # Return Value:
    ans
}


# ------------------------------------------------------------------------------




gfuncSlider <- 
    function(B = 10)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Slider for generator function, density and probability

    # FUNCTION:

    # Graphic Frame:
    par(mfrow = c(2, 2), cex = 0.7)

    # Internal Function:
    refresh.code = function(...)
    {
        # Startup Counter:
        .counter <- setRmetricsOptions(".counter") + 1
        setRmetricsOptions(.counter = .counter)
        if (.counter < 6) return ()

        # Sliders:
        Copula = as.integer(.sliderMenu(no = 1))
        type = ellipticalList()
        type = type[Copula]
        Type = c("Normal", "Cauchy", "Student-t", "Logistic", "Laplace",
            "Kotz", "Exponential Power")
        Type = Type[Copula]
        N = .sliderMenu(no = 2)
        nu = .sliderMenu(no = 3)
        r = .sliderMenu(no = 4)
        s = .sliderMenu(no = 5)
        rho = .sliderMenu(no = 6)
        L = 6.5

        # Parameters:
        param = NULL
        if (Copula == 3) param = nu
        if (Copula == 6) param = r
        if (Copula == 7) param = c(r, s)
        prefactor = gfunc(param = param, type = type)[[1]]
        Lambda = as.character(round(prefactor, digits = 3))
        Nu = R = S = NA
        if (Copula == 3) Nu = as.character(round(nu, digits = 1))
        if (Copula >= 6) R = as.character(round(r, digits = 1))
        if (Copula == 7) S = as.character(round(s, digits = 1))
        delta = 10/N

        # Bivariate Density:
        x = y = seq(-4, 4, length = 101)
        D = delliptical2d(grid2d(x), rho = rho, param = param,
            type = type, output = "list")

        # Plot 1:
        Limit = ""
        if (Copula == 3 & nu == 1) Limit = "| [Cauchy]"
        if (Copula == 6 & r == 1) Limit = "| [Normal]"
        if (Copula == 7 & s == 1) Limit = "| [Kotz]"
        if (Copula == 7 & r == 1 & s == 1) Limit = "| [Normal]"
        lambda = gfunc(param = param, type = type)
        x = seq(0, L, length = N)
        y = gfunc(x, param = param, type = type)
        y.norm = gfunc(x, type = "norm")
        plot(x, y, type = "l", ylab = "g", ylim = c(0, 1))
        abline(h = 0, lty = 3, col = "grey")
        lines(x, y.norm, lty = 3, col = "red")
        title(main = paste("Generator:", Type, Limit, "\nPre-Factor:", Lambda))
        mtext("Dotted Curve: Normal Generator", side = 4, col = "grey",
            cex = 0.7)

        # Plot 2 - Density:
        x = seq(-L, L, length = N)
        y = .delliptical(x, param = param, type = type)
        y.norm = .delliptical(x, type = "norm")
        plot(x, y, type = "l", ylab = "Density", ylim = c(0, 0.65))
        abline(h = 0, lty = 3, col = "grey")
        abline(v = 0, lty = 3, col = "grey")
        lines(x, y.norm, lty = 3, col = "red")
        Y = 2*integrate(.delliptical, 0, Inf, param = param, type = type)[[1]]
        Y = as.character(round(Y, 2))
        .velliptical = function(x, param, type) x^2*.delliptical(x, param, type)
        V = 2*integrate(.delliptical, 0, Inf, param = param, type = type)[[1]]
        V = as.character(round(V, 2))
        mtext(paste("Normalization Test:", Y, " |  Variance Test:", V),
            side = 4, col = "grey", cex = 0.7)
        if (type == "t") {
            title(main = paste(Type, "Density\n nu =", Nu))
        } else if (type == "kotz") {
            title(main = paste(Type, "Density\n r =", R))
        } else if (type == "epower") {
            title(main = paste(Type, "Density\n r =", R, "s =", S))
        } else {
            title(main = paste(Type, "Density\n "))
        }

        # Plot 3 - Probability:
        x = seq(-L, L, length = N)
        y = .pelliptical(x, param = param, type = type)
        y.norm = .pelliptical(x, type = "norm")
        plot(x, y, type = "l", ylab = "Probability", ylim = c(0, 1))
        abline(h = 0, lty = 3, col = "grey")
        abline(h = 1, lty = 3, col = "grey")
        abline(h = 0.5, lty = 3, col = "grey")
        lines(x, y.norm, lty = 3, col = "red")
        p95 = .qelliptical(0.95, param = param, type = type)
        P95 = as.character(round(p95, digits = 2))
        abline(v = p95, lty = 3)
        abline(v = -p95, lty = 3)
        q95 = .pelliptical(p95, param = param, type = type)
        points(+p95, q95, pch = 19, cex = 0.5)
        points(-p95, 1-q95, pch = 19, cex = 0.5)
        mtext("Dots: Probability(Quantile(0.95)) Test", side = 4,
            col = "grey", cex = 0.7)
        Title = paste(Type, "Probability\n 95% =", P95)
        title(main = Title)

        # Plot 4 - Bivariate Density:
        contour(D, levels = c(0.001, 0.01, 0.025, 0.05, 0.1),
            xlab = "x", ylab = "y")
        title(main = paste("Bivariate Density\nrho = ", as.character(rho)))
        grid()

        # Reset Frame:
        par(mfrow = c(2, 2), cex = 0.7)
    }

    # Open Slider Menu:
    setRmetricsOptions(.counter = 0)
    .sliderMenu(refresh.code,
        names       = c("Copula",   "N", "3: nu",  "6|7: r",  "7: s", "rho"),
        minima      = c(       1,    50,       1,       0.1,     0.1, -0.95),
        maxima      = c(       7,  2000,       B,         B,       B,  0.95),
        resolutions = c(       1,    50,     0.1,       0.1,     0.1,  0.05),
        starts      = c(       1,   100,       4,         1,       1,  0.00))
}


# ------------------------------------------------------------------------------


.pelliptical <- 
    function(q, param = NULL, type = ellipticalList(),
    alternative = TRUE, subdivisions = 100)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Probability function for univariate elliptical distributions

    # Arguments:
    #   q -  a numeric vector
    #   param - NULL, a numeric value, or a numeric vector adding.
    #       additional parameters to the generator function.
    #   type -  a character string denoting the type of distribution.
    #       This may be either
    #       "norm" for the normal distribution, or
    #       "cauchy" for the Cauchy distribution, or
    #       "t" for the Student-t distribution, or
    #       "logistic" for the logistic distribution, or
    #       "laplace" for the distribution, or
    #       "kotz" for the original Kotz distribution, or
    #       "epower" for the exponential power distribution.

    # Details:
    #   The probability is computed by integration using the generator
    #       function. If an alternative faster algorithm is available,
    #       this one is used by default.

    # FUNCTION:

    # Type:
    type = match.arg(type)

    # Alternative Available?
    if (type == "logistic") alternative = FALSE
    if (type == "laplace") alternative = FALSE
    if (type == "kotz") alternative = FALSE
    if (type == "epower") alternative = FALSE

    # Original Function:
    # Fq1 = function (x, Q, param, type) {
    #    acos(abs(Q)/sqrt(x)) * gfunc(x, param, type) }
    # Transformed Function: u = exp(-x+Q^2)
    Fq2 =
    function (x, Q, param, type)
    {
        Q^2 * acos(sqrt(x))/x^2 * gfunc(Q^2/x, param, type)
    }

    # Add Default Parameters:
    if (is.null(param)) {
        if (type == "t") param = c(nu = 4)
        if (type == "kotz") param = c(r = 1)
        if (type == "epower") param = c(r = 1, s = 1)
    }

    # Probability:
    ans = NULL
    if (alternative) {
        ans = NA
        if (type[1] == "norm") ans = pnorm(q)
        if (type[1] == "cauchy") ans = pt(q, df = 1) # pcauchy(q)
        if (type[1] == "t") ans = pt(q, df = param[[1]])
        if (type[1] == "kotz") ans = dnorm(q, sd = 1/sqrt(param[[1]]))
    } else {
        lambda = gfunc(param = param, type = type)[[1]]
        ans = NULL
        for ( Q in q ) {
            # More Precise Adaptive Rule:
            # p = lambda * integrate(Fq1, lower = Q^2, upper = Inf, Q = Q,
            # param = param, type = type, subdivisions = subdivisions)[[1]]
            p = lambda*integrate(Fq2, lower = .Machine$double.eps^0.5,
                upper = 1, Q = Q, param = param, type = type,
                stop.on.error = FALSE, subdivisions = subdivisions)[[1]]
            if (Q > 0) p = 1 - p
            if (abs(Q) < .Machine$double.eps^0.5) p = 0.5
            ans = c(ans, p)
        }
    }

    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


.delliptical <- 
    function(x, param = NULL, type = ellipticalList(), alternative = TRUE,
    subdivisions = 100)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Density function for univariate elliptical distributions

    # Arguments:
    #   x -  a numeric vector
    #   param - NULL, a numeric value, or a numeric vector adding.
    #       additional parameters to the generator function.
    #   type -  a character string denoting the type of distribution.
    #       This may be either
    #       "norm" for the normal distribution, or
    #       "cauchy" for the Cauchy distribution, or
    #       "t" for the Student-t distribution, or
    #       "logistic" for the logistic distribution, or
    #       "laplace" for the distribution, or
    #       "kotz" for the original Kotz distribution, or
    #       "epower" for the exponential power distribution.
    #   alternative -  a logical flag. Should alternatively used a
    #       faster algorithm if available? By default TRUE.

    # Details:
    #   The density is computed by integration using the generator
    #       function. If an alternative faster algorithm is available,
    #       this one is used by default.

    # FUNCTION:

    # Type:
    type = match.arg(type)

    # Alternative Available?
    if (type == "logistic") alternative = FALSE
    if (type == "laplace") alternative = FALSE
    if (type == "kotz") alternative = FALSE
    if (type == "epower") alternative = FALSE

    # Original Function:
    # fq1 = function (x, Q, param, type) {
    #    gfunc(x, param, type) / ( sqrt(x - Q^2) ) }
    # Transformed Function: log(x)^2 = x - Q^2
    fq2 = function (x, Q, param, type) {
        2 * gfunc(log(x)^2+Q^2, param, type) / x  }


    # Add Default Parameters:
    if (is.null(param)) {
        if (type == "t") param = c(nu = 4)
        if (type == "kotz") param = c(r = 1)
        if (type == "epower") param = c(r = 1, s = 1)
    }

    # Normalizing constant lambda:
    lambda = gfunc(param = param, type = type)[[1]]

    # Density:
    ans = NULL
    if (alternative) {
        ans = NA
        if (type[1] == "norm") ans = dnorm(x)
        if (type[1] == "cauchy") ans = dt(x, df = 1) # dcauchy(x)
        if (type[1] == "t") ans = dt(x, df = param[[1]])
        if (type[1] == "kotz") ans = dnorm(x, sd = 1/sqrt(param[[1]]))
    } else {
        lambda = gfunc(param = param, type = type)[[1]]
        ans = NULL
        for ( Q in x ) {
            # More Precise Adaptive Rule:
            # p = lambda*integrate(fq1, lower = Q^2, upper = Inf, Q = Q,
            #     param = param, type = type)[[1]]
            p = lambda*integrate(fq2, lower = 0, upper = 1, Q = Q, param =
                param, type = type, stop.on.error = FALSE,
                subdivisions = subdivisions)[[1]]
            ans = c(ans, p)
        }
    }

    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


.qelliptical <- 
    function(p, param = NULL, type = ellipticalList(), alternative = TRUE)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Quantile function for univariate elliptical distributions

    # Arguments:
    #   p -  a numeric vector
    #   param - NULL, a numeric value, or a numeric vector adding.
    #       additional parameters to the generator function.
    #   type -  a character string denoting the type of distribution.
    #       This may be either
    #       "norm" for the normal distribution, or
    #       "cauchy" for the Cauchy distribution, or
    #       "t" for the Student-t distribution, or
    #       "logistic" for the logistic distribution, or
    #       "laplace" for the distribution, or
    #       "kotz" for the original Kotz distribution, or
    #       "epower" for the exponential power distribution.
    #   alternative - a logical flag. Should be an alternative
    #       faster algorithm used and not the standard algorithm?

    # Details:
    #   The probability is computed by integration using the generator
    #       function. If an alternative faster algorithm is available,
    #       this one is used by default.

    # FUNCTION:

    # Type:
    type = match.arg(type)

    # Alternative Available?
    if (type == "laplace") alternative = FALSE
    if (type == "kotz") alternative = FALSE
    if (type == "epower") alternative = FALSE

    # Add Default Parameters:
    if (is.null(param)) {
        if (type == "t") param = c(nu = 4)
        if (type == "kotz") param = c(r = 1)
        if (type == "epower") param = c(r = 1, s = 1)
    }

    # Probability:
    ans = NULL
    if (alternative) {
        ans = NA
        if (type[1] == "norm") ans = qnorm(p)
        if (type[1] == "cauchy") ans = qcauchy(p)
        if (type[1] == "t") ans = qt(p, df = param[[1]])
        if (type[1] == "logistic") ans = .qlogistic(p)
        if (type[1] == "kotz") ans = dnorm(p, sd = 1/sqrt(param[[1]]))
    } else {
        froot <-
        function(x, p, param, type)
        {
            .pelliptical(q = x, param = param, type = type) - p
        }
        ans = NULL
        for (pp in p) {
            if (pp < .Machine$double.eps) {
                ans = c(ans, -Inf)
            } else if (pp > 1 - .Machine$double.eps) {
                ans = c(ans, Inf)
            } else {
                lower = -1
                upper = +1
                counter = 0
                iteration = NA
                while (is.na(iteration)) {
                    iteration = .unirootNA(f = froot, interval = c(lower,
                        upper), param = param, type = type, p = pp)
                    counter = counter + 1
                    lower = lower - 2^counter
                    upper = upper + 2^counter
                }
                ans = c(ans, iteration)
            }
        }
    }

    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


.qlogistic <- 
    function(p)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Fast Quantile function for the logistic distribution

    # FUNCTION:

    # Table:
    data = .qlogisticTable

    # Quantiles:
    P = (sign(p-1/2)+1)/2 - sign(p-1/2)*p
    ans = sign(0.5-p) * approx(x = data[, 2], y = data[, 1], xout = P)$y

    # p Boundary:
    index = which(p < 0.001 & p > 0)
    if (length(index) > 0) {
        ans[index] =
            .qelliptical(p[index], type = "logistic", alternative = FALSE)
    }
    index = which(p > 1-0.001 & p < 1)
    if (length(index) > 0) {
        ans[index] =
            .qelliptical(p[index], type = "logistic", alternative = FALSE)
    }
    ans[p == 0.5] = 0
    ans[p == 0] = -Inf
    ans[p == 1] = Inf

    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


.qlogisticData <- 
    function (dump = FALSE )
{   
    # A function implemented by Diethelm Wuertz

    # FUNCTION:

    # Range:
    p = seq(0.001, 0.500, by = 0.001)

    # Quantiles by Integration:
    froot =
    function(x, p)
    {
        .pelliptical(x, type = "logistic") - p
    }

    X = NULL
    for (P in p) {
        lower = -1
        upper = +1
        counter = 0
        iteration = NA
        while (is.na(iteration)) {
            iteration = .unirootNA(f = froot, interval = c(lower, upper), p = P)
            counter = counter + 1
            lower = lower - 2^counter
            upper = upper + 2^counter
        }
        X = c(X, iteration)
    }
    Y = .pelliptical(X, type = "logistic")
    .qlogisticTable = data.frame(cbind(X = X, Y = Y))

    # Dump:
    if (dump) dump(".qlogisticTable", "qlogisticTable.R")

    # Return Value:
    invisible(.qlogisticTable)
}


# ------------------------------------------------------------------------------


.qlogisticTable <- 
    structure(list(
    X = c(
    -3.28961095698868, -3.08838952917050, -2.96495324742154,
    -2.87441959067841, -2.80235793855428, -2.74216585623189, -2.69027685632636,
    -2.64454429353653, -2.60362855984489, -2.56644066721983, -2.53234562858188,
    -2.50082506289166, -2.47148229772502, -2.44400686160488, -2.41815107184679,
    -2.39371408491880, -2.37053072103299, -2.34846344243532, -2.32739647429067,
    -2.30723141378254, -2.28791218251300, -2.26930432832626, -2.25137887299253,
    -2.23407971956745, -2.21735785753918, -2.20116988800091, -2.18547719458243,
    -2.17024525936846, -2.15544310253510, -2.14104280951074, -2.12701913726668,
    -2.11334918152059, -2.1000120943421,  -2.08698884326631, -2.07426200479765,
    -2.06181558657013, -2.04963487351506, -2.03770629424047, -2.02601730451810,
    -2.01455628530110, -2.00331245315227, -1.99227557959128, -1.98143672734491,
    -1.97078698069461, -1.96031819468915, -1.95002274531123, -1.93989348557452,
    -1.92992370611904, -1.92010709998667, -1.91043773071894, -1.90091000365456,
    -1.89151863995147, -1.88225865304357, -1.87312532726254, -1.86411419838938,
    -1.85522103592924, -1.84644182692950, -1.83777276118156, -1.82921021766652,
    -1.82074871277337, -1.81238487394136, -1.80411802330877, -1.79594515630526,
    -1.78786340687312, -1.77987003898311, -1.77196243800415, -1.76413810662673,
    -1.75639465334193, -1.74872978943490, -1.7412023575414,  -1.73368818564377,
    -1.7262462929639,  -1.71887474483312, -1.71157168326167, -1.70433532304102,
    -1.69716394791141, -1.69005590701564, -1.68300961157971, -1.67602353180290,
    -1.66909619394163, -1.66222617757269, -1.65541211302264, -1.64865267895148,
    -1.64194660007940, -1.63529264504669, -1.62868962439725, -1.62213638867753,
    -1.61563182664270, -1.60917486356300, -1.60276445962361, -1.59639952641875,
    -1.59007925325713, -1.58380261455573, -1.57756869582143, -1.57137660985193,
    -1.56522549889593, -1.55911452857102, -1.55304289016109, -1.54700979879211,
    -1.54101449249394, -1.53505623130807, -1.52915261112904, -1.52326458881694,
    -1.51741167285209, -1.5115931890518,  -1.50580848264567, -1.50005691749585,
    -1.49433787536073, -1.48865075519602, -1.48299497249102, -1.47736995863785,
    -1.47177516036924, -1.46621003903524, -1.4606740702944,  -1.4551667434371,
    -1.44968756093365, -1.4442360379646,  -1.43881170197258, -1.43341409223487,
    -1.42804275945499, -1.42269726537273, -1.41737718239130, -1.41208209321679,
    -1.40681159053369, -1.40156527665306, -1.39634276321712, -1.39114367089538,
    -1.38596762909750, -1.38081427569799, -1.37568325677215, -1.37057422634254,
    -1.36548684613562, -1.36042078534799, -1.35537572042174, -1.35035133482848,
    -1.34534731886170, -1.34036336943694, -1.33539918989953, -1.33045448983943,
    -1.32552898491301, -1.32062239667119, -1.31573445239396, -1.31086488493075,
    -1.30601343254650, -1.30117983877317, -1.29636385226647, -1.29156522666752,
    -1.28678372046926, -1.28201909688745, -1.27727112373605, -1.27253957330672,
    -1.26782422225068, -1.26312485147296, -1.25844124601302, -1.25377319494646,
    -1.24912049128124, -1.24448293185896, -1.23986031725926, -1.23525245170728,
    -1.23065914298397, -1.22608020233923, -1.22151544440770, -1.21696225226631,
    -1.21242200176641, -1.20789558040343, -1.20338280567676, -1.19888349860929,
    -1.19439748365017, -1.18992458855333, -1.18546464439883, -1.18101748533607,
    -1.17658294861366, -1.17216087446918, -1.16775110605060, -1.16341452449662,
    -1.15902890823835, -1.15465514386589, -1.15029304960735, -1.14594255361642,
    -1.14160347928037, -1.13727568813146, -1.13295904407152, -1.12865341331407,
    -1.12435866432809, -1.12007466778373, -1.11580129649967, -1.11153842461814,
    -1.10728593065946, -1.10304369280653, -1.09881159197813, -1.09458951100203,
    -1.09037733457144, -1.08617494920272, -1.08198224319439, -1.07779910658730,
    -1.07362543112605, -1.06946111022151, -1.06528196623425, -1.06113353941502,
    -1.05699432071886, -1.05286419954172, -1.04874306739911, -1.04463081785862,
    -1.04058838163140, -1.0364935858852,  -1.03240736512243, -1.02832962049820,
    -1.02426025492196, -1.02019917300514, -1.01614628101095, -1.01210148680617,
    -1.00806469981498, -1.00403583097456, -1.00001479269255, -0.995987692143912,
    -0.99198188126047, -0.987982882923126,-0.983990646036694,-0.980005119840238,
    -0.976026253905185,-0.972053998133449,-0.968088302755545,-0.964129118328723,
    -0.96023743089134, -0.956291121335992,-0.952351176345156,-0.948417547764324,
    -0.944491447079358,-0.940578633091236,-0.936672395329904,-0.932772667777382,
    -0.928879385227229,-0.924992483271505,-0.921111898287973,-0.917237567427541,
    -0.913369428601936,-0.909507420471605,-0.905651482433843,-0.901801554611138,
    -0.89795757783973, -0.894119493658381,-0.890287244297353,-0.886460772667583,
    -0.882640022350059,-0.878824937585387,-0.875015463263552,-0.871211544913856,
    -0.867413128695051,-0.863620161385636,-0.859832590374336,-0.856050363650754,
    -0.852273429796184,-0.848501737974597,-0.844735237923773,-0.840973879946612,
    -0.837217614902577,-0.833466394199305,-0.82972016978436, -0.825978894137128,
    -0.82224252026086, -0.81851100167485, -0.81478429240676, -0.81106234698506,
    -0.807345120431618,-0.803632568254409, -0.799924646440352, -0.796221311448274,
    -0.79252252020199, -0.788828230083503, -0.785138398926326, -0.781452985008911,
    -0.777771947048196,-0.774095244193254, -0.770422836019059, -0.766754682520353,
    -0.763090744605994,-0.759430982092061, -0.755775356696648, -0.752123830034272,
    -0.748476364110065,-0.744832921314231, -0.741193464432255, -0.737557956577102,
    -0.733926361277322,-0.730298642409806, -0.726674764210157, -0.723054691267668,
    -0.719438388520393,-0.715825821250292, -0.712216955078455, -0.708611755960407,
    -0.705010190181483,-0.701412224352282, -0.697817825404193, -0.69422696058499,
    -0.690639597454501,-0.687055703880342, -0.68347524803372, -0.679898198385306,
    -0.676324523701167,-0.672754193038763, -0.669187175743011, -0.665623441442406,
    -0.662062960045204,-0.658505701735664, -0.654951636970346, -0.651400736474472,
    -0.647852971238333,-0.644308312513762, -0.640766731810652, -0.637228200893533,
    -0.633692691778194,-0.630160176728364, -0.626630628252438, -0.623104019100251,
    -0.619580322259904,-0.616059510954638, -0.612541558639745, -0.60902643899954,
    -0.60551412594436, -0.602004593607625, -0.598497816342926, -0.594993768721168,
    -0.591492425527744,-0.587993761759759, -0.584497752623287, -0.581004373530672,
    -0.577513600097867,-0.574025408141805, -0.570539773677817, -0.567056672917079,
    -0.563576082264099,-0.560097978314237, -0.55662233785126,  -0.553149137844933,
    -0.549678355448641,-0.546209967997043, -0.542743953003765, -0.539280288159113,
    -0.53581895132783, -0.532359920546873, -0.528903174023226, -0.525448690131743,
    -0.521996447413011,-0.518546424571258, -0.515098600472271, -0.511652954141354,
    -0.508209464761305,-0.504768111670429, -0.501328874360566, -0.497891732475152,
    -0.494456665807302,-0.491023654297919, -0.487592678033829, -0.484163717245933,
    -0.480736752307392,-0.477311763731828, -0.473888732171551, -0.470467638415803,
    -0.467048463389035,-0.463631188149193, -0.460215793886032, -0.45680226191945,
    -0.453390573697839,-0.449980710796464, -0.44657265491585,  -0.443166387880198,
    -0.43976189163582, -0.436359148249581, -0.432958139907375, -0.429558848912609,
    -0.426161257684707,-0.422765348757635, -0.419371104778434, -0.415978508505784,
    -0.41258754280857, -0.409198190664474, -0.405810435158577, -0.402424259481984,
    -0.399039646930456,-0.395656580903062, -0.392275044900847, -0.388895022525509,
    -0.385516497478098,-0.382139453557724, -0.378763874660278, -0.375389744777172,
    -0.372017047994090,-0.368645768489747, -0.36527589053467,  -0.361907398489987,
    -0.358540276806227,-0.355174510022136, -0.351810082763504, -0.348446979742001,
    -0.345085185754029,-0.341724685679586, -0.338365464481134, -0.335007507202486,
    -0.331650798967702,-0.328295324979995, -0.324941070520646, -0.321588020947929,
    -0.318236161696055,-0.314885478274112, -0.311535956265027, -0.308187581324529,
    -0.304840339180130,-0.301494152790354, -0.298149133027608, -0.294805203648722,
    -0.291462350657759,-0.28812056012555,  -0.284779814605342, -0.281440107419829,
    -0.278101421295547,-0.274763742560770, -0.271427057605776, -0.268091352881914,
    -0.264756614900680,-0.261422830232804, -0.258089985507338, -0.254758067410761,
    -0.251427062686079,-0.248096958131943, -0.244767740601768, -0.241439397002857,
    -0.238111914295537,-0.234785279492298, -0.231459479656942, -0.228134501903726,
    -0.224810333396537,-0.221486961348039, -0.218164373018853, -0.214842555716734,
    -0.21152149679575, -0.208201183655463, -0.204881603740137, -0.201562744537922,
    -0.198244593580064,-0.194927138440111, -0.191610366733128, -0.188294266114916,
    -0.184978824281230,-0.181664028967013, -0.178349867945625, -0.175036329028080,
    -0.171723400062291,-0.168411068932309, -0.165099323557576, -0.161788151892181,
    -0.158477541924117,-0.15516748167454,  -0.151857959197035, -0.148548962576890,
    -0.145240479930364,-0.141932499403962, -0.138625009173723, -0.135317997444490,
    -0.132011452449203,-0.128705362448187, -0.125399715728434, -0.122094500602899,
    -0.118789705409793,-0.115485318511868, -0.112181328295712, -0.108877723171037,
    -0.105574491569969,  -0.102271621946344,  -0.0989691027749957,
    -0.0956669225510829, -0.0923650697894468, -0.0890635330240574,
    -0.0857623008076265, -0.0824613617115197, -0.0791607043262014,
    -0.0758603172625963, -0.0725601891549516, -0.069260308666142,
    -0.0659606644967098, -0.0626612453993505, -0.0593620402006167,
    -0.0560630378306748, -0.0527642273584164, -0.0494655980205171,
    -0.0461671392155751, -0.0428688404062201, -0.0395706908379783,
    -0.0362726789608858, -0.0329747914456716, -0.0296770116896015,
    -0.0263793175678826, -0.0230816777192755, -0.0197840449696766,
    -0.0164872991738574, -0.0132014769537292, -0.0099153529555198,
    -0.00656807273015348,-0.00328316565264927, 0),
    Y = c(0.000999989790249075,
    0.00199997659631066,0.0030000533744429, 0.00400032302535487,
    0.00500028183011033,0.00600024440622486, 0.00700036400791461,
    0.00800057651215374,0.00899926075702241, 0.00999936918645765,
    0.0109994715833541,0.0119995614818643,0.0129996374478306,0.0139997002917028,
    0.0149997517033347,0.0159997935471271,0.0169998275539801,0.0179998552300887,
    0.0189998778128955,0.0199998963093886,0.0209984200243254,0.0219986195324741,
    0.0229987985092580,0.0239989556602773,0.0249990929774641,0.0259992125331238,
    0.0269993163154686,0.0279994063127642,0.028999484230876, 0.0299995516446606,
    0.0309996099536415,0.0319996603894939,0.0329997040283535,0.033999741805388,
    0.034999774530007, 0.0359998029007155,0.0369998275190071,0.0379998489022468,
    0.0389998674948231,0.0399998836784007,0.0409998977807025,0.0419999100844623,
    0.0429999208323128,0.0439999302262356,0.0449999384488791,0.0459999456542143,
    0.0469999519731455,0.0479999575239014,0.0489999624033934,0.0499999666971179,
    0.050999970479087, 0.0519999738134476,0.0529999767558638,0.0539999793547017,
    0.0549999816520199,0.0559999836844612,0.0569999854839726,0.0579999870784435,
    0.0589999884922518,0.0600002322660306,0.0610007382842464,0.0620012176556511,
    0.0630016721028811,0.0640021032232844,0.0650025124989912,0.0660029014066227,
    0.0670032710226466,0.0680036226385918,0.0690039573610388,0.0699961909442808,
    0.0709964157890431,0.0719966272564782,0.0729968261635939,0.0739970132882051,
    0.074997189351493, 0.075997355027219, 0.076997510945194, 0.0779976576944622,
    0.0789977958262262,0.0799979258565364,0.0809980482687644,0.0819981635158824,
    0.0829982720225615,0.0839983741871086,0.0849984703832535,0.0859985609618,
    0.0869986462521522,0.0879987265637271,0.0889988021872626,0.0899988733960295,
    0.0909989402880648,0.0919990034241728,0.0929990628713735,0.0939991188435381,
    0.0949991716757273,0.0959992212888924,0.0969992679961586,0.0979993119655162,
    0.0989993533551862,0.0999993923142102,0.100999428983003, 0.101996361340610,
    0.102996667458178, 0.103996948294911, 0.104997205893694, 0.105997442136406,
    0.106997658756622, 0.107997857350875, 0.108998039389096, 0.109998206224314,
    0.110998359101641, 0.111998499159902, 0.112998627466711, 0.113998744983746,
    0.114998852601866, 0.115998951139948, 0.116999041350541, 0.117999123925120,
    0.118999199498883, 0.119999268655243, 0.120999331929934, 0.121999389814826,
    0.122999442762193, 0.123999491184829, 0.124999535463806, 0.125999575947988,
    0.126999612957384, 0.127999646785496, 0.12899967770146,  0.129999705952036,
    0.13099973176343,  0.131999755342970, 0.132999776880658, 0.133999796550584,
    0.134999814512239, 0.135999830911716, 0.136999845882814, 0.137999859548057,
    0.138999872019624, 0.139999883400211, 0.140999893783817, 0.141999903256468,
    0.142999911896882, 0.143999919777082, 0.144999926962952, 0.145999933514754,
    0.146999939487602, 0.147999944931889, 0.148999949893692, 0.149999954415129,
    0.150999958534699, 0.151999962287943, 0.152999965706274, 0.153999968819446,
    0.154999971654293, 0.155999974235324, 0.156999976584922, 0.157999978723527,
    0.158999980669799, 0.159999982440775, 0.160999984052003, 0.161999985517674,
    0.163000522714305, 0.164001257335553, 0.165001955709628, 0.166002619635127,
    0.167003250824463, 0.168003850907805, 0.169004421443053, 0.170004963894409,
    0.171005479673739, 0.17200597011836,  0.173006436500641, 0.174006880031024,
    0.174993398219765, 0.175993761862996, 0.176994106194154, 0.177994440428791,
    0.178994740610310, 0.179995032746045, 0.180995309240802, 0.181995570901728,
    0.182995818495872, 0.183996052752091, 0.184996274362880, 0.185996483986112,
    0.186996682428530, 0.187996869918613, 0.188997047203267, 0.189997214817756,
    0.190997373270305, 0.191997523043414, 0.192997664595105, 0.193997798360121,
    0.194997924751063, 0.195998044159472, 0.196998156956866, 0.198004064106562,
    0.199004781577514, 0.200005456363768, 0.201006090989097, 0.202006687834809,
    0.203007249147298, 0.203992879823315, 0.204993345345794, 0.205993781542165,
    0.20699419018612,  0.207994572949000, 0.208994931405369, 0.209995267038304,
    0.210995581244424, 0.211995875338658, 0.212996150558767, 0.213996408069636,
    0.215000093310879, 0.216000369275620, 0.217000821976959, 0.218001445611404,
    0.219002234487624, 0.220003183024853, 0.221004285751305, 0.222005537302608,
    0.223006932420246, 0.223992987554550, 0.22499462733946,  0.225996395713614,
    0.226998287825741, 0.227999977350380, 0.228999972825903, 0.229999967831223,
    0.230999962364020, 0.231999956424305, 0.232999950014237, 0.233999943137955,
    0.234999935801416, 0.235999928012238, 0.236999919779558, 0.237999911113890,
    0.238999902026997, 0.23999989253177,  0.240999882642106, 0.241999872372804,
    0.242999861739457, 0.243999850758354, 0.244999839446394, 0.245999827820992,
    0.246999815900002, 0.247999803701641, 0.248999791244413, 0.249999778547048,
    0.250999765628433, 0.251999752507560, 0.252999739203466, 0.253999725735183,
    0.254999712121693, 0.255999698381881, 0.256999684534495, 0.257999670598112,
    0.258999656591096, 0.259999642531572, 0.260999628437393, 0.261999614326116,
    0.262999600214975, 0.263999586120861, 0.264999572060296, 0.265999558049426,
    0.266999544103992, 0.267999530239325, 0.268999516470329, 0.269999502811470,
    0.270999489276768, 0.271999475879788, 0.272999462633633, 0.273999449550939,
    0.274999436643870, 0.275999423928137, 0.276999411406967, 0.277999399095057,
    0.278999387002665, 0.279999375139568, 0.280999363515069, 0.281999352133693,
    0.282999341012350, 0.283999330154682, 0.284999319568118, 0.285999309259627,
    0.286999299235725, 0.287999289502482, 0.288999280065525, 0.289999270930049,
    0.290999262100821, 0.291999253582187, 0.292999245378082, 0.293999237492035,
    0.294999229927180, 0.295999222686262, 0.296999215771647, 0.297999209185332,
    0.298999202928950, 0.299999197003782, 0.300999191410765, 0.301999186150503,
    0.302999181223271, 0.303999176629033, 0.304999172367443, 0.305999168437859,
    0.306999164839350, 0.307999161570707, 0.308999158630451, 0.309999156016846,
    0.3109991537279,   0.311999151761384, 0.312999150114833, 0.313999148785561,
    0.314999147770665, 0.315999147067038, 0.316999146671377, 0.317999146580187,
    0.318999146789796, 0.319999147296359, 0.320999148095870, 0.321999149184165,
    0.322999150556936, 0.323999152209734, 0.324999154137979, 0.325999156336969,
    0.326999158801883, 0.327999161527793, 0.328999164509669, 0.329999167742386,
    0.330999171220733, 0.331999174939416, 0.332999178893067, 0.333999183076251,
    0.334999187483469, 0.335999192109171, 0.336999196947752, 0.337999201993567,
    0.338999207240932, 0.33999921268413,  0.340999218317417, 0.341999224135028,
    0.342999230131179, 0.343999236300076, 0.344999242635916, 0.345999249132897,
    0.346999255785214, 0.347999262587073, 0.348999269532686, 0.349999276616284,
    0.350999283832112, 0.351999291174441, 0.352999298637566, 0.353999306215812,
    0.354999313903537, 0.355999321695133, 0.356999329585035, 0.357999337567717,
    0.3589993456377,   0.359999353789552, 0.360999362017892, 0.361999370317391,
    0.362999378682777, 0.363999387108836, 0.364999395590413, 0.365999404122415,
    0.366999412699814, 0.367999421317649, 0.368999429971024, 0.369999438655116,
    0.37099944736517,  0.371999456096506, 0.372999464844517, 0.373999473604672,
    0.374999482372517, 0.375999491143676, 0.376999499913850, 0.377999508678821,
    0.378999517434455, 0.379999526176695, 0.380999534901569, 0.381999543605189,
    0.382999552283749, 0.383999560933529, 0.384999569550893, 0.38599957813229,
    0.386999586674258, 0.387999595173416, 0.388999603626474, 0.389999612030226,
    0.390999620381554, 0.391999628677426, 0.392999636914899, 0.393999645091113,
    0.394999653203298, 0.395999661248771, 0.396999669224934, 0.397999677129276,
    0.398999684959373, 0.399999692712886, 0.400999700387563, 0.401999707981236,
    0.402999715491823, 0.403999722917325, 0.40499973025583,  0.405999737505508,
    0.406999744664612, 0.407999751731479, 0.408999757699432, 0.409999764564645,
    0.410999771332861, 0.411999778002731, 0.412999784572988, 0.413999792115251,
    0.414999798496769, 0.415999804775648, 0.416999810950938, 0.417999817021771,
    0.41899982298735,  0.419999828846956, 0.420999834599940, 0.421999840245729,
    0.422999845783821, 0.423999851213786, 0.424999856535265, 0.425999861747966,
    0.42699986685167,  0.427999871846223, 0.428999876731538, 0.429999881507595,
    0.43099988617444,  0.431999890732181, 0.432999895180991, 0.433999899521106,
    0.434999903752824, 0.435999907876499, 0.436999911892552, 0.437999915801457,
    0.438999919603749, 0.439999923300017, 0.440999926890908, 0.441999930377125,
    0.442999933759421, 0.443999937038606, 0.444999940215541, 0.445999943291138,
    0.446999946266358, 0.447999949142211, 0.448999951919758, 0.449999954600106,
    0.450999957184407, 0.45199995967386,  0.452999962069707, 0.453999964373235,
    0.454999966585772, 0.455999968708689, 0.456999970743396, 0.457999972691344,
    0.458999974554023, 0.45999997633296,  0.460999978029718, 0.461999979645896,
    0.462999981183131, 0.46399998264309,  0.464999984027476, 0.465999985338022,
    0.466999986576494, 0.467999987744687, 0.468999988844425, 0.469999989877562,
    0.470999990845976, 0.471999991751573, 0.472999992596279, 0.473999993382044,
    0.474999994110834, 0.47599999478463,  0.476999995405422, 0.477999995975202,
    0.478999996495962, 0.479999996969683, 0.480999997398343, 0.481999997783931,
    0.48299999812849,  0.483999998434203, 0.484999998703553, 0.485999998939509,
    0.486999999145696, 0.487999999326349, 0.488999999485856, 0.489999999628004,
    0.490999999755448, 0.49199999986999,  0.492999999973083, 0.494000000064158,
    0.494999711104453, 0.495996057848887, 0.496992418350526, 0.49800732088262,
    0.499003719695633, 0.499999999819624)),
    .Names = c("X", "Y"),
    row.names = c("1",
      "2",  "3",  "4",  "5",  "6",  "7", "8",  "9",  "10",  "11", "12", "13",
     "14", "15", "16", "17", "18", "19", "20", "21", "22",  "23", "24",
     "25", "26", "27", "28", "29", "30", "31", "32", "33",  "34", "35",
     "36", "37", "38", "39", "40", "41", "42", "43", "44",  "45", "46",
     "47", "48", "49", "50", "51", "52", "53", "54", "55",  "56", "57",
     "58", "59", "60", "61", "62", "63", "64", "65", "66",  "67", "68",
     "69", "70", "71", "72", "73", "74", "75", "76", "77",  "78", "79",
     "80", "81", "82", "83", "84", "85", "86", "87", "88",  "89", "90",
     "91", "92", "93", "94", "95", "96", "97", "98", "99",  "100",
    "101", "102", "103", "104", "105", "106", "107", "108", "109",
    "110", "111", "112", "113", "114", "115", "116", "117", "118",
    "119", "120", "121", "122", "123", "124", "125", "126", "127",
    "128", "129", "130", "131", "132", "133", "134", "135", "136",
    "137", "138", "139", "140", "141", "142", "143", "144", "145",
    "146", "147", "148", "149", "150", "151", "152", "153", "154",
    "155", "156", "157", "158", "159", "160", "161", "162", "163",
    "164", "165", "166", "167", "168", "169", "170", "171", "172",
    "173", "174", "175", "176", "177", "178", "179", "180", "181",
    "182", "183", "184", "185", "186", "187", "188", "189", "190",
    "191", "192", "193", "194", "195", "196", "197", "198", "199",
    "200", "201", "202", "203", "204", "205", "206", "207", "208",
    "209", "210", "211", "212", "213", "214", "215", "216", "217",
    "218", "219", "220", "221", "222", "223", "224", "225", "226",
    "227", "228", "229", "230", "231", "232", "233", "234", "235",
    "236", "237", "238", "239", "240", "241", "242", "243", "244",
    "245", "246", "247", "248", "249", "250", "251", "252", "253",
    "254", "255", "256", "257", "258", "259", "260", "261", "262",
    "263", "264", "265", "266", "267", "268", "269", "270", "271",
    "272", "273", "274", "275", "276", "277", "278", "279", "280",
    "281", "282", "283", "284", "285", "286", "287", "288", "289",
    "290", "291", "292", "293", "294", "295", "296", "297", "298",
    "299", "300", "301", "302", "303", "304", "305", "306", "307",
    "308", "309", "310", "311", "312", "313", "314", "315", "316",
    "317", "318", "319", "320", "321", "322", "323", "324", "325",
    "326", "327", "328", "329", "330", "331", "332", "333", "334",
    "335", "336", "337", "338", "339", "340", "341", "342", "343",
    "344", "345", "346", "347", "348", "349", "350", "351", "352",
    "353", "354", "355", "356", "357", "358", "359", "360", "361",
    "362", "363", "364", "365", "366", "367", "368", "369", "370",
    "371", "372", "373", "374", "375", "376", "377", "378", "379",
    "380", "381", "382", "383", "384", "385", "386", "387", "388",
    "389", "390", "391", "392", "393", "394", "395", "396", "397",
    "398", "399", "400", "401", "402", "403", "404", "405", "406",
    "407", "408", "409", "410", "411", "412", "413", "414", "415",
    "416", "417", "418", "419", "420", "421", "422", "423", "424",
    "425", "426", "427", "428", "429", "430", "431", "432", "433",
    "434", "435", "436", "437", "438", "439", "440", "441", "442",
    "443", "444", "445", "446", "447", "448", "449", "450", "451",
    "452", "453", "454", "455", "456", "457", "458", "459", "460",
    "461", "462", "463", "464", "465", "466", "467", "468", "469",
    "470", "471", "472", "473", "474", "475", "476", "477", "478",
    "479", "480", "481", "482", "483", "484", "485", "486", "487",
    "488", "489", "490", "491", "492", "493", "494", "495", "496",
    "497", "498", "499", "500"),
    class = "data.frame"
)


################################################################################

