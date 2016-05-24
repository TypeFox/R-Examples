
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
# FUNCTION:                 EXTREME VALUE COPULAE PARAMETER:
#  evList                    Returns list of implemented extreme value copulae
#  evParam                   Sets Default parameters for an extreme value copula
#  evCheck                   Checks if parameters are in the valid range
#  evRange                   Returns the range of valid parameter values
# FUNCTION:                 EXTREME VALUE COPULAE GENERATOR FUNCTION:
#  Afunc                     Computes Dependence function
#  AfuncSlider               Displays interactively dependence function
#  .AfuncFirstDer             Computes Derivative of dependence function
#  .AfuncSecondDer            Computes 2nd Derivative of dependence function
################################################################################


################################################################################
# FUNCTION:                 EXTREME VALUE COPULAE PARAMETER:
#  evList                    Returns list of implemented extreme value copulae
#  evParam                   Sets parameters for an extreme value copula
#  evRange                   Returns the range of valid parameter values
#  evCheck                   Checks if parameters are in the valid range


evList =
function()
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Returns list of implemented extreme value copulae

    # Compose List:
    ans = c("gumbel", "galambos", "husler.reiss", "tawn", "bb5")

    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


evParam =
function(type = evList())
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Sets default parameters for extreme value copulae

    # Arguments:
    #   type - a character string naming the copula. By default the
    #       "gumbel" copula will be chosen.

    # Value:
    #   returns a list with two elements, 'param' sets the parameters
    #       which may be a vector, 'range' the range with minimum and
    #       maximum values for each of the parameters. For the "pi" and
    #       "m" copula NULL will be returned.

    # FUNCTION:

    # Settings:
    type = match.arg(type)
    ans = list(copula = type)

    # Select:
    if ( type == "gumbel" ) {
        ans$param = c(delta = 2)
        ans$range = c(1, Inf) }
    if ( type == "galambos" ) {
        ans$param = c(delta = 2)
        ans$range = c(0, Inf) }
    if ( type == "husler.reiss" ) {
        ans$param = c(delta = 2)
        ans$range = c(0, Inf) }
    if ( type == "tawn" ) {
        ans$param = c(alpha = 2, beta = 1/2, r = 2)
        ans$range = c(0, 1, 0, 1, 1, Inf) }
    if ( type == "bb5" ) {
        ans$param = c(delta = 2, theta = 2)
        ans$range = c(0, Inf, 0, Inf) }

    # Some more, yet untested and undocumented:
    if ( type == "gumbelII" ) {
        ans$param = c(alpha = 2)
        ans$range = NULL }
    if ( type == "marshall.olkin" ) {
        ans$param = c(alpha1 = 2, alpha2 = 2)
        ans$range = NULL }
    if ( type == "pi" ) {
        ans$param = NULL
        ans$range = NULL }
    if ( type == "m" ) {
        ans$param = NULL
        ans$range = NULL }

    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


evRange =
function(type = evList())
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Returns the range of valid parameter values

    # Examples:
    #   evRange("galambos")
    #   evRange("bb5")

    # FUNCTION:

    # Type:
    type = match.arg(type)

    # Range:
    ans = evParam(type)$range
    Names1 = rep(c("lower", "upper"), times = length(ans)/2)
    Names2 = rep(names(evParam(type)$param), each = 2)
    names(ans) = paste(Names1, Names2, sep = ".")
    attr(ans, "control")<-type

    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


evCheck =
function(param, type = evList())
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Checks if parameters are in the valid range

    # FUNCTION:

    # Type:
    type = match.arg(type)

    # Check
    range = evRange(type)
    nParam = length(range)/2
    j = -1
    J = 0
    for (i in 1:nParam) {
        j = j + 2
        J = J + 2
        if (param[i] < range[j] | param[i] > range[J]) {
            print(c(param = param[i]))
            print(c(range = c(range[j], range[J])))
            stop("param is out of range")
        }
    }

    # Return Value:
    invisible(TRUE)
}


################################################################################
# FUNCTION:                  EXTREME VALUE COPULAE GENERATOR FUNCTION:
#  Afunc                      Computes Dependence function
#  AfuncSlider                Displays interactively dependence function
#  .AfuncFirstDer              Computes Derivative of dependence function
#  .AfuncSecondDer             Computes 2nd Derivative of dependence function


Afunc =
function(x, param = NULL, type = evList())
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Computes dependence function for extreme value copulae

    # Arguments:
    #   x - a numeric vector, with values ranging between
    #       zero and one
    #   param - numeric parameter vector, if set to NULL then
    #       default values are taken
    #   type - character string naming the type of copula,
    #       by default "gumbel"

    # Details:
    #   Extreme Value Copulae can be represented in the form
    #
    #      C(u,v) = exp { log(uv)*A[log(u)/log(uv)] }
    #
    #   where A:[0,1] -> [1/2,1] is a convex function
    #   such that max(x,1-x) < A(x) < 1 for all x in [0,1].

    # Notes:
    #   Copulae included also in EVANESCE:
    #       gumbel, galambos, husler.reiss, tawn, bb5
    #   Additionally - not yet tested and documented
    #       gumbelII, marshall.olkin, pi[Cperp], m[Cplus]

    # References:
    #   Bouye E. (2000), Copulas for Finance: A Reading Guide and
    #       Some Applications, (see the Table on page 49).
    #   Insightful Corp, EVANESCE Implementation in S-PLUS
    #       FinMetrics Module.

    # FUNCTION:

    # Missing x:
    if (missing(x)) x = (0:10)/10

    # Type:
    type = type[1]
    if (is.null(param)) param = evParam(type)$param
    names(param) = names(evParam(type)$param)

    # Compute Dependence Function:
    if (type == "gumbel") {
        # 1 <= alpha < Inf
        alpha = param[1]
        if (alpha == 1) A = rep(1, times = length(x)) else
        A = (x^alpha + (1-x)^alpha)^(1/alpha)
    }
    if (type == "galambos") {
        # 0 <= alpha < Inf
        alpha = param[1]
        A = 1 - (x^(-alpha) + (1-x)^(-alpha))^(-1/alpha)
    }
    if (type == "husler.reiss") {
        # 0 <= alpha <= Inf
        alpha = param[1]
        A = x * pnorm(1/alpha + 0.5*alpha*log(x/(1-x))) +
            (1-x) * pnorm(1/alpha - 0.5*alpha*log(x/(1-x)))
    }
    if (type == "tawn") {
        # 0 <= alpha <=1
        # 0 <= beta <= 1
        # 1 <= r < Inf
        alpha = param[1]
        beta = param[2]
        r = param[3]
        if (alpha == 0 | beta == 0 | r == 1) A = rep(1, times = length(x)) else
        A = 1 - beta +(beta-alpha)*x + ( (alpha*x)^r + (beta*(1-x))^r )^(1/r)
    }
    if (type == "bb5") {
        # 0 < delta < Inf
        # 1 <= theta Inf
        delta = param[1]
        theta = param[2]
        if (theta == 1) return(Afunc(x, param, "galambos")) else
        A = ( x^theta + (1-x)^theta -
            ( x^(-delta*theta) + (1-x)^(-delta*theta) )^(-1/delta))^(1/theta)
    }

    # Some more, yet untested and undocumented:
    if (type == "gumbelII") {
        # 0 <= alpha < Inf
        alpha = param[1]
        A = alpha*x^2 - alpha*x + 1
    }
    if (type == "marshall.olkin") {
        alpha1 = param[1]
        alpha2 = param[2]
        A = NULL
        for (i in 1:length(x)) A = c(A, max(1-alpha1*x[i], 1-alpha2*(1-x[i])))
    }
    if (type == "pi" || type == "Cperp") {
        # No parameters
        A = rep(1, times = length(x))
    }
    if (type == "m" || type == "Cplus") {
        # No parameters
        A = NULL
        for (i in 1:length(x)) A = c(A, max(x[i], 1-x[i]))
    }

    # Result:
    attr(A, "control") <- unlist(list(param = param, type = type))

    # Return Value:
    A
}


# ------------------------------------------------------------------------------


AfuncSlider =
function()
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Displays interactively the dependence function

    # Graphic Frame:
    par(mfrow = c(2, 2), cex = 0.7)

    # Internal Function:
    refresh.code = function(...)
    {
        # Startup Counter:
        .counter <- getRmetricsOptions(".counter") + 1
        setRmetricsOptions(.counter = .counter)
        if (.counter < 10) return ()

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

        # Plot A:
        plot(x = (0:N)/N, Afunc(x = (0:N)/N, param = param, type = type),
            ylim = c(0.5, 1), type = "l", xlab = "x", ylab = "A", main = Title)
        lines(c(0.0, 1.0), c(1.0, 1.0), col = "steelblue", lty = 3)
        lines(c(0.0, 0.5), c(1.0, 0.5), col = "steelblue", lty = 3)
        lines(c(0.5, 1.0), c(0.5, 1.0), col = "steelblue", lty = 3)
        points(x = c(0, 1), Afunc(x = c(0, 1), param = param, type = type),
            col = "red")
        # Plot A':
        plot(x = (0:N)/N, .AfuncFirstDer(x = (0:N)/N, param = param, type = type),
            type = "l", xlab = "x", ylab = "A'", main = Title)
        points(x = c(0, 1),
            .AfuncFirstDer(x = c(0, 1), param = param, type = type), col = "red")
        # Plot A'':
        plot(x = (0:N)/N, .AfuncSecondDer(x = (0:N)/N, param = param, type = type),
            type = "l", xlab = "x", ylab = "A''", main = Title)
        points(x = c(0, 1),
            .AfuncSecondDer(x = c(0, 1), param = param, type = type), col = "red")

        # Reset Frame:
        par(mfrow = c(2, 2), cex = 0.7)
    }

    # Open Slider Menu:
    setRmetricsOptions(.counter = 0)
    C = c("Gumbel: delta", "Galambos: delta", "Husler-Reis: delta",
          "Tawn: alpha", "... beta", "... r", "BB5: delta", "... theta")
    .sliderMenu(refresh.code,
        names =       c("Copula", "N", C), #gal hr  tawn               bb5
        minima =      c(1,   100,  1.0, 0.00, 0.00, 0.00, 0.00, 1.0, 0.0, 1.0),
        maxima =      c(5, 10000, 10.0, 10.0, 10.0, 1.00, 1.00, 10., 10., 10.),
        resolutions = c(1,   100, 0.05, 0.05, 0.05, 0.01, 0.01, 0.1, 0.1, 0.1),
        starts =      c(1,  5000, 1.00, 0.00, 0.00, 0.00, 0.00, 1.0, 0.0, 1.0))
}


# ------------------------------------------------------------------------------


.AfuncFirstDer =
function(x, param = NULL, type = evList(), eps = 1.0e-6 )
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   #   Computes derivaive of dependence function

    # Arguments:
    #   x - a numeric vector, with values ranging between
    #       zero and one
    #   param - numeric parameter vector, if set to NULL then
    #       default values are taken
    #   type - character string naming the type of copula,
    #       by default "gumbel"

    # Details:
    #   Extreme Value Copulae can be represented in the form
    #
    #      C(u,v) = exp { log(uv)*A[log(u)/log(uv)] }
    #
    #   where A:[0,1] -> [1/2,1] is a convex function
    #   such that max(x,1-x) < A(x) < 1 for all x in [0,1].

    # Notes:
    #   Copulae included also in EVANESCE:
    #       gumbel, galambos, husler.reiss, tawn, bb5
    #   Additionally - not yet tested and documented
    #       gumbelII, marshall.olkin, pi[Cperp], m[Cplus]

    # References:
    #   Bouye E. (2000), Copulas for Finance: A Reading Guide and
    #       Some Applications, (see the Table on page 49).
    #   Insightful Corp, EVANESCE Implementation in S-PLUS
    #       FinMetrics Module.

    # FUNCTION:

    # Missing x:
    if (missing(x)) x = (0:10)/10

    # Type:
    type = type[1]
    if (is.null(param)) param = evParam(type)$param
    names(param) = names(evParam(type)$param)

    # Settings for Maple Output:
    Pi = pi
    ln = function(x) { log(x) }
    erf = function (x) { 2*pnorm(sqrt(2)*x)-1 }

    # Compute Derivative:
    if (type == "gumbel") {
        # alpha >= 1
        alpha = param[1]
        # Maple Generated Output:
        if (alpha == 1) A1 = rep(0, times = length(x)) else {
        A1 =
        (x^alpha+(1-x)^alpha)^(1/alpha)/alpha*(x^alpha*alpha/x-(1-x)^alpha*
        alpha/(1-x))/(x^alpha+(1-x)^alpha)
        A1[x < eps] = -1
        A1[x > 1-eps] = 1 }
    }
    if (type == "galambos") {
        # 0 <= alpha < Inf
        alpha = param[1]
        # Maple Generated Output:
        if (alpha == 0) A1 = rep(1, times = length(x)) else {
        A1 =
        (x^(-alpha)+(1-x)^(-alpha))^(-1/alpha)/alpha*(-x^(-alpha)*alpha/x+(
        1-x)^(-alpha)*alpha/(1-x))/(x^(-alpha)+(1-x)^(-alpha))
        A1[x < eps] = -1
        A1[x > 1-eps] = 1 }
    }
    if (type == "husler.reiss") {
        # 0 <= alpha <= Inf
        alpha = param[1]
        # Maple Generated Output:
        if (alpha == 0) A1 = rep(1, times = length(x)) else {
        A1 =
        .5*erf(1/2*(1/alpha+.5*alpha*ln(x/(1-x)))*2^(1/2))+.2500000000/Pi^(
        1/2)*exp(-1/2*(1/alpha+.5*alpha*ln(x/(1-x)))^2)*alpha*(1/(1-x)+x/(1
        -x)^2)*(1-x)*2^(1/2)-.5*erf(1/2*(1/alpha-.5*alpha*ln(x/(1-x)))*2^(1
        /2))-.2500000000*(1-x)^2/Pi^(1/2)*exp(-1/2*(1/alpha-.5*alpha*ln(x/(
        1-x)))^2)*alpha*(1/(1-x)+x/(1-x)^2)/x*2^(1/2)
        A1[x < eps] = -1
        A1[x > 1-eps] = 1 }
    }
    if (type == "tawn") {
        # 0 <= alpha < Inf
        # beta <= 1
        # 1 <= r < Inf
        alpha = param[1]
        beta = param[2]
        r = param[3]
        # Maple Generated Output:
        if (alpha == 0 | beta == 0 | r == 1) A1 = rep(0, length(x)) else {
        A1 =
        beta-alpha+((alpha*x)^r+(beta*(1-x))^r)^(1/r)/r*((alpha*x)^r*r/x-(
        beta*(1-x))^r*r/(1-x))/((alpha*x)^r+(beta*(1-x))^r)
        A1[x < eps] = -alpha
        A1[x > 1-eps] = beta }
    }
    if (type == "bb5") {
        # 0 < delta < Inf
        # 1 <= theta < Inf
        delta = param[1]
        theta = param[2]
        # Maple Generated Output:
        if (theta == 1) return(.AfuncFirstDer(x, param, "galambos")) else
        A1 = (x^theta+(1-x)^theta-(x^(-delta*theta)+(1-x)^(-delta*theta))^(-1/
            delta))^(1/theta)/theta*(x^theta*theta/x-(1-x)^theta*theta/(1-x)+(x
            ^(-delta*theta)+(1-x)^(-delta*theta))^(-1/delta)/delta*(-x^(-delta*
            theta)*delta*theta/x+(1-x)^(-delta*theta)*delta*theta/(1-x))/(x^(-
            delta*theta)+(1-x)^(-delta*theta)))/(x^theta+(1-x)^theta-(x^(-delta
            *theta)+(1-x)^(-delta*theta))^(-1/delta))
        A1[x < eps] = -1
        A1[x > 1-eps] = 1
    }

    # Some more, yet untested and undocumented:
    if (type == "gumbelII") {
        # 0 <= alpha < Inf
        alpha = param[1]
        A1 = 2*alpha*x-alpha
    }
    if (type == "marshall.olkin") {
        alpha1 = param[1]
        alpha2 = param[2]
        A1 = NULL
        for (i in 1:length(x)) {
            if (x[i] < 0) A1 = c(A1, -alpha1)
            if (x[i] > 0) A1 = c(A1, alpha2)
            if (x[i] == 0) A1 = c(A1, NA) }
    }
    if (type == "pi" || type == "Cperp") {
        A1 = rep(0, times = length(x))
    }
    if (type == "m" || type == "Cplus") {
        A1 = sign(x-1/2)
    }

    # Result:
    attr(A1, "control") <- unlist(list(param = param, type = type))

    # Return Value:
    A1
}


# ------------------------------------------------------------------------------


.AfuncSecondDer =
function(x, param = NULL, type = evList())
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Computes 2nd derivative of dependence function

    # Arguments:
    #   x - a numeric vector, with values ranging between
    #       zero and one
    #   param - numeric parameter vector, if set to NULL then
    #       default values are taken
    #   type - character string naming the type of copula,
    #       by default "gumbel"

    # Details:
    #   Extreme Value Copulae can be represented in the form
    #
    #      C(u,v) = exp { log(uv)*A[log(u)/log(uv)] }
    #
    #   where A:[0,1] -> [1/2,1] is a convex function
    #   such that max(x,1-x) < A(x) < 1 for all x in [0,1].

    # Note:
    #   The five Copulae considered in EVANESCE are:
    #       gumbel, galambos, husler.reis, tawn, bb5
    #   Furthermore, added are:
    #       pi|Cperp, gumbelII, marshall.olkin, m|Cplus

    # References:
    #   Bouye E. (2000), Copulas for Finance: A Reading Guide and
    #       Some Applications, (see the Table on page 49).
    #   Insightful Corp, EVANESCE Implementation in S-PLUS
    #       FinMetrics Module.

    # FUNCTION:

    # Missing x:
    if (missing(x)) x = (0:10)/10

    # Type:
    type = type[1]
    if (is.null(param)) param = evParam(type)$param
    names(param) = names(evParam(type)$param)

    # Settings for Maple Output:
    Pi = pi
    ln = function(x) { log(x) }
    erf = function (x) { 2*pnorm(sqrt(2)*x)-1 }

    # Compute 2nd Derivative:
    if (type == "gumbel") {
        # alpha >= 1
        alpha = param[1]
        # Maple Generated Output:
        if (alpha == 1) A2 = rep(0, times = length(x)) else
        A2 = (x^alpha+(1-x)^alpha)^(1/alpha)/alpha^2*(x^alpha*alpha/x-(1-x)^
            alpha*alpha/(1-x))^2/(x^alpha+(1-x)^alpha)^2+(x^alpha+(1-x)^alpha)^
            (1/alpha)/alpha*(x^alpha*alpha^2/x^2-x^alpha*alpha/x^2+(1-x)^alpha*
            alpha^2/(1-x)^2-(1-x)^alpha*alpha/(1-x)^2)/(x^alpha+(1-x)^alpha)-(x
            ^alpha+(1-x)^alpha)^(1/alpha)/alpha*(x^alpha*alpha/x-(1-x)^alpha*
            alpha/(1-x))^2/(x^alpha+(1-x)^alpha)^2
    }
    if (type == "galambos") {
        # 0 <= alpha < Inf
        alpha = param[1]
        # Maple Generated Output:
        if (alpha == 0) A2 = rep(0, times = length(x)) else
        if (alpha == 1) A2 = rep(2, times = length(x)) else
        A2 = -(x^(-alpha)+(1-x)^(-alpha))^(-1/alpha)/alpha^2*(-x^(-alpha)*alpha/
            x+(1-x)^(-alpha)*alpha/(1-x))^2/(x^(-alpha)+(1-x)^(-alpha))^2+(x^(-
            alpha)+(1-x)^(-alpha))^(-1/alpha)/alpha*(x^(-alpha)*alpha^2/x^2+x^(
            -alpha)*alpha/x^2+(1-x)^(-alpha)*alpha^2/(1-x)^2+(1-x)^(-alpha)*
            alpha/(1-x)^2)/(x^(-alpha)+(1-x)^(-alpha))-(x^(-alpha)+(1-x)^(-
            alpha))^(-1/alpha)/alpha*(-x^(-alpha)*alpha/x+(1-x)^(-alpha)*alpha/
            (1-x))^2/(x^(-alpha)+(1-x)^(-alpha))^2
    }
    if (type == "husler.reiss") {
        # 0 <= alpha <= Inf
        alpha = param[1]
        # Maple Generated Output:
        if (alpha == 0) A2 = rep(0, times = length(x)) else
        A2 = .2500000000/Pi^(1/2)*exp(-1/2*(1/alpha+.5*alpha*ln(x/(1-x)))^2)*
            alpha*(1/(1-x)+x/(1-x)^2)/x*(1-x)*2^(1/2)-.1250000000/Pi^(1/2)*(1/
            alpha+.5*alpha*ln(x/(1-x)))*alpha^2*(1/(1-x)+x/(1-x)^2)^2/x*(1-x)^2*
            exp(-1/2*(1/alpha+.5*alpha*ln(x/(1-x)))^2)*2^(1/2)+.2500000000/Pi^(
            1/2)*exp(-1/2*(1/alpha+.5*alpha*ln(x/(1-x)))^2)*alpha*(2/(1-x)^2+2
            *x/(1-x)^3)*(1-x)*2^(1/2)-.2500000000/Pi^(1/2)*exp(-1/2*(1/alpha+.5
            *alpha*ln(x/(1-x)))^2)*alpha*(1/(1-x)+x/(1-x)^2)*2^(1/2)+.75000000/
            Pi^(1/2)*exp(-1/2*(1/alpha-.5*alpha*ln(x/(1-x)))^2)*alpha*(1/(1-x
            )+x/(1-x)^2)/x*(1-x)*2^(1/2)-.1250000000*(1-x)^3/Pi^(1/2)*(1/alpha-
            .5*alpha*ln(x/(1-x)))*alpha^2*(1/(1-x)+x/(1-x)^2)^2/x^2*exp(-1/2*(1
            /alpha-.5*alpha*ln(x/(1-x)))^2)*2^(1/2)-.2500000000*(1-x)^2/Pi^(1/2
            )*exp(-1/2*(1/alpha-.5*alpha*ln(x/(1-x)))^2)*alpha*(2/(1-x)^2+2*x/(
            1-x)^3)/x*2^(1/2)+.2500000000*(1-x)^2/Pi^(1/2)*exp(-1/2*(1/alpha-.5*
            alpha*ln(x/(1-x)))^2)*alpha*(1/(1-x)+x/(1-x)^2)/x^2*2^(1/2)
    }
    if (type == "tawn") {
        # 0 <= alpha, beta <= 1, 1 <= r < Inf
        alpha = param[1]
        beta = param[2]
        r = param[3]
        # Maple Generated Output:
        if (alpha == 0 | beta == 0 | r == 1) A2 = rep(0, length(x)) else
        A2 = ((alpha*x)^r+(beta*(1-x))^r)^(1/r)/r^2*((alpha*x)^r*r/x-(beta*(1-x)
            )^r*r/(1-x))^2/((alpha*x)^r+(beta*(1-x))^r)^2+((alpha*x)^r+(beta*(1
            -x))^r)^(1/r)/r*((alpha*x)^r*r^2/x^2-(alpha*x)^r*r/x^2+(beta*(1-x))
            ^r*r^2/(1-x)^2-(beta*(1-x))^r*r/(1-x)^2)/((alpha*x)^r+(beta*(1-x))^
            r)-((alpha*x)^r+(beta*(1-x))^r)^(1/r)/r*((alpha*x)^r*r/x-(beta*(1-x
            ))^r*r/(1-x))^2/((alpha*x)^r+(beta*(1-x))^r)^2
        # A2[x<1e-12] = 0
        # A2[x>1-1e-12] = 0
    }
    if (type == "bb5") {
        # delta > 0, theta >= 1
        delta = param[1]
        theta = param[2]
        # Maple Generated Output:
        if (theta == 1) return(.AfuncSecondDer(x, param, "galambos")) else
        A2 = (x^theta+(1-x)^theta-(x^(-delta*theta)+(1-x)^(-delta*theta))^(-1/
            delta))^(1/theta)/theta^2*(x^theta*theta/x-(1-x)^theta*theta/(1-x)+
            (x^(-delta*theta)+(1-x)^(-delta*theta))^(-1/delta)/delta*(-x^(-
            delta*theta)*delta*theta/x+(1-x)^(-delta*theta)*delta*theta/(1-x))/
            (x^(-delta*theta)+(1-x)^(-delta*theta)))^2/(x^theta+(1-x)^theta-(x^
            (-delta*theta)+(1-x)^(-delta*theta))^(-1/delta))^2+(x^theta+(1-x)^
            theta-(x^(-delta*theta)+(1-x)^(-delta*theta))^(-1/delta))^(1/theta)/
            theta*(x^theta*theta^2/x^2-x^theta*theta/x^2+(1-x)^theta*theta^2/(
            1-x)^2-(1-x)^theta*theta/(1-x)^2-(x^(-delta*theta)+(1-x)^(-delta*
            theta))^(-1/delta)/delta^2*(-x^(-delta*theta)*delta*theta/x+(1-x)^(
            -delta*theta)*delta*theta/(1-x))^2/(x^(-delta*theta)+(1-x)^(-delta*
            theta))^2+(x^(-delta*theta)+(1-x)^(-delta*theta))^(-1/delta)/delta*
            (x^(-delta*theta)*delta^2*theta^2/x^2+x^(-delta*theta)*delta*theta/
            x^2+(1-x)^(-delta*theta)*delta^2*theta^2/(1-x)^2+(1-x)^(-delta*
            theta)*delta*theta/(1-x)^2)/(x^(-delta*theta)+(1-x)^(-delta*theta))
            -(x^(-delta*theta)+(1-x)^(-delta*theta))^(-1/delta)/delta*(-x^(-
            delta*theta)*delta*theta/x+(1-x)^(-delta*theta)*delta*theta/(1-x))^
            2/(x^(-delta*theta)+(1-x)^(-delta*theta))^2)/(x^theta+(1-x)^theta-(
            x^(-delta*theta)+(1-x)^(-delta*theta))^(-1/delta))-(x^theta+(1-x)^
            theta-(x^(-delta*theta)+(1-x)^(-delta*theta))^(-1/delta))^(1/theta)/
            theta*(x^theta*theta/x-(1-x)^theta*theta/(1-x)+(x^(-delta*theta)+(
            1-x)^(-delta*theta))^(-1/delta)/delta*(-x^(-delta*theta)*delta*
            theta/x+(1-x)^(-delta*theta)*delta*theta/(1-x))/(x^(-delta*theta)+(
            1-x)^(-delta*theta)))^2/(x^theta+(1-x)^theta-(x^(-delta*theta)+(1-x
            )^(-delta*theta))^(-1/delta))^2
    }

    # Some more, yet untested and undocumented:
    if (type == "gumbelII") {
        alpha = param[1]
        A2 = rep(2*alpha, times = length(x))
    }
    if (type == "marshall.olkin") {
        alpha1 = param[1]
        alpha2 = param[2]
        A2 = rep(0, times = length(x))
    }
        if (type == "pi" || type == "Cperp") {
        A2 = rep(0, times = length(x))
    }
    if (type == "m" || type == "Cplus") {
        A2 = rep(0, times = length(x))
    }

    # Result:
    attr(A2, "control") <- unlist(list(param = param, type = type))

    # Return Value:
    A2
}


################################################################################

