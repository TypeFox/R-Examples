
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
# FUNCTION                   KENDALL'S TAU AND SPEARMAN'S RHO:
#  evTau                      Returns Kendall's tau for extreme value copulae
#  .ev1Tau                     Computes Kendall's tau from dependency function
#  .ev2Tau                     Computes Kendall's tau from integration
#  evRho                      Returns Spearman's rho for extreme value copulae
#  .ev1Rho                     Computes Spearman's rho from dependency function
#  .ev2Rho                     Computes Spearman's rho from integration
# FUNCTION:                  EXTREME VALUE COPULAE TAIL DEPENDENCE:
#  evTailCoeff                Computes tail dependence for extreme value copulae
#  evTailCoeffSlider          Plots extreme value tail dependence function
################################################################################


################################################################################
# FUNCTION                   KENDALL'S TAU AND SPEARMAN'S RHO:
#  evTau                      Returns Kendall's tau for extreme value copulae
#  evRho                      Returns Spearman's rho for extreme value copulae


evTau =
function(param = NULL, type = evList(), alternative = FALSE)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Computes Kendall's tau for an extreme value copula

    # Example:
    #   evTau(alternative = FALSE)
    #   evTau(alternative = TRUE)

    # FUNCTION:

    # Kendall's Tau:
    if (!alternative) {
        # Default Method:
        ans = .ev1Tau(param, type)
    } else {
        # Alternative Method:
        ans = .ev2Tau(param, type)
    }

    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


.ev1Tau =
function(param = NULL, type = evList())
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Computes Kendall's tau from dependency function

    # FUNCTION:

    # Type:
    type = match.arg(type)

    # Default Parameters:
    if (is.null(param)) param = evParam(type)$param

    # Kendall's Tau Integrand:
    fun = function(x, param, type) {
        # To be integrated from 0 to 1 ...
        A = Afunc(x = x, param = param, type = type)
        A2 = .AfuncSecondDer(x, param, type)
        f = (x*(1-x)/A) * A2
        f
    }

    # Get control attribute from:
    attribute = Afunc(0.5, param, type)

    # Integrate:
    ans = integrate(fun, 0, 1, param = param, type = type)
    Tau = c(Tau = ans[[1]])

    # Add Control Attribute:
    attr(Tau, "control")<-attr(attribute, "control")

    # Return Value:
    Tau
}


# ------------------------------------------------------------------------------


.ev2Tau =
function(param = NULL, type = evList())
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Computes Kendall's tau from integration

    # Example:
    #   .ev2Tau()

    # FUNCTION:

    # Type:
    type = match.arg(type)

    # Default Parameters:
    if (is.null(param)) param = evParam(type)$param

    # Kendall's Tau Minus Rho/3 Double Integrand:
    fun = function(x, y, ...) {
        D = devCopula(x, y, alternative = FALSE, ...)
        D[is.na(D)] = 0
        f = 4 *
            ( pevCopula(x, y, alternative = FALSE, ...) - x*y) * D
        f
    }

    # Get control attribute from:
    attribute = Afunc(0.5, param, type)

    # Integrate:
    ans = integrate2d(fun, param = param, type = type, error = 1e-8)
    Tau = c(Tau = ans[[1]] + .ev2Rho(param, type)/3)

    # Add Control Attribute:
    attr(Tau, "control")<-attr(attribute, "control")

    # Return Value:
    Tau
}


# ------------------------------------------------------------------------------


evRho =
function(param = NULL, type = evList(), alternative = FALSE)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Computes Spearman's rho for an extreme value copula

    # Example:
    #   evRho(alternative = FALSE)
    #   evRho(alternative = TRUE)

    # FUNCTION:

    # Spearman's Rho:
    if (!alternative) {
        # Default Method:
        ans = .ev1Rho(param, type)
    } else {
        # Alternative Method:
        ans = .ev2Rho(param, type)
    }

    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


.ev1Rho =
function(param = NULL, type = evList())
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Computes Spearman's rho from dependency function

    # Example:
    #   .ev1Rho()

    # FUNCTION:

    # Type:
    type = match.arg(type)

    # Default Parameters:
    if (is.null(param)) param = evParam(type)$param

    # Spearman's Rho Integrand:
    fun = function(x, param, type) {
        # To be integrated from 0 to 1 ...
        A = Afunc(x = x, param = param, type = type)
        f = ( 12 / (A+1)^2 ) - 3
        f
    }

    # Get control attribute from:
    attribute = Afunc(0.5, param, type)

    # Integrate:
    ans = integrate(fun, 0, 1, param = param, type = type)
    Rho = c(Rho = ans[[1]])

    # Add Control Attribute:
    attr(Rho, "control")<-attr(attribute, "control")

    # Return Value:
    Rho
}


# ------------------------------------------------------------------------------


.ev2Rho =
function(param = NULL, type = evList())
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Computes Spearman's rho from integration

    # Example:
    #   .ev2Rho()

    # FUNCTION:

    # Type:
    type = match.arg(type)

    # Default Parameters:
    if (is.null(param)) param = evParam(type)$param

    # Spearman's Rho Integrand:
    fun = function(x, y, ...) {
        f = 12 * pevCopula(x, y, ...) - 3
        f
    }

    # Get control attribute from:
    attribute = Afunc(0.5, param, type)

    # Integrate:
    ans = integrate2d(fun, param = param, type = type)
    Rho = c(Rho = ans[[1]])

    # Add Control Attribute:
    attr(Rho, "control")<-attr(attribute, "control")

    # Return Value:
    Rho
}


################################################################################
# FUNCTION:                  EXTREME VALUE COPULAE TAIL DEPENDENCE:
#  evTailCoeff                Computes tail dependence for extreme value copulae
#  evTailCoeffSlider          Plots extreme value tail dependence function


evTailCoeff =
function(param = NULL, type = evList())
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Tail Dependence for extreme value copulae

    # Example:
    #   evTailCoeff()

    # FUNCTION:

    # Type:
    type = match.arg(type)

    # Default Parameters:
    if (is.null(param)) param = evParam(type)$param

    # Limit:
    lambdaU = 2-2*Afunc(0.5, param, type)[[1]]
    lambdaL = 0
    ans = c(lower = lambdaL, upper = lambdaU)

    # Add Control Attribute:
    attr(ans, "control") <-
        unlist(list(copula = "ev", param = param, type = type))

    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


evTailCoeffSlider =
function(B = 10)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Displays interactively perspective plots of tail coefficient

    # Example:
    #   evTailCoeffSlider()

    # FUNCTION:

    # Graphic Frame:
    par(mfrow = c(1, 1))

    # Internal Function:
    refresh.code = function(...)
    {
        # Startup:
        .counter <- getRmetricsOptions(".counter") + 1
        setRmetricsOptions(.counter = .counter)
        if (.counter < 10) return()

        # Sliders:
        Type = evList()
        Copula = .sliderMenu(no = 1)
        N = .sliderMenu(no = 2)
        if (Copula <= 3)
            param = c(delta = .sliderMenu(no = Copula + 2))
        if (Copula == 4)
            param = c(alpha = .sliderMenu(no = 6),
                beta = .sliderMenu(no = 7), r = .sliderMenu(no = 8))
        if (Copula == 5)
            param = c(delta = .sliderMenu(no = 9), theta = .sliderMenu(no = 10))

        # Title:
        type = Type[Copula]
        subTitle = paste(paste(names(param) , "="), param, collapse = " | " )
        Title = paste(" ", type, "\n", subTitle)

        # Plot:
        u = seq(0, 0.5, length = N+1)[-1]
        C.uu = pevCopula(u, u, param, type)
        lambda = C.uu/u
        v = seq(0.5, 1, length = N+1)[-(N+1)]
        C.uu = pevCopula(v, v, param, type)
        lambda = c(lambda, (1-2*v+C.uu)/(1-v))
        x = c(u, v)
        plot(x, lambda, xlim = c(0, 1), ylim = c(0, 1),
            pch = 19, col = "steelblue", xlab = "u")
        title(main = Title)
        grid()

        # Add Points:
        points(x = 0, y = 0, pch = 19, col = "red")
        points(x = 1, y = 2-2*Afunc(0.5, param, type), pch = 19, col = "red")

        # Lines:
        abline(h = 0, col = "grey")
        abline(v = 0.5, col = "grey")

        # Reset Frame:
        par(mfrow = c(1, 1))
    }

    setRmetricsOptions(.counter = 0)
    # Open Slider Menu:
    C = c("1 Gumbel: delta", "2 Galambos: delta", "3 Husler-Reis: delta",
          "4 Tawn: alpha", "... beta", "... r", "5 BB5: delta", "... theta")
    .sliderMenu(refresh.code,
        names = c("Copula", "N", C),
                         #   N gumbel galamb h.r  tawn-tawn-tawn  bb5-bb5
        minima =      c(1,  10,     1,    0,   0,    0,   0,   1,   0,  1),
        maxima =      c(5, 100,     B,    B,   B,    1,   1,   B,   B,  B),
        resolutions = c(1,  10,   .05,  .05, .05,  .01, .01,  .1,  .1, .1),
        starts =      c(1,  20,     2,    1,   1,   .5,  .5,   2,   1,  2))
}


################################################################################

