
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
# FUNCTION:                  EXTREME VALUE COPULAE RANDOM VARIATES:
#  revCopula                  Generates extreme value copula random variates
#  revSlider                  isplays interactively plots of random variates
# FUNCTION:                  EXTREME VALUE COPULAE PROBABILIY:
#  pevCopula                  Computes extreme value copula probability
#  pevSlider                  Displays interactively plots of probability
#  .pev1Copula                 EV copula probability via dependence function
#  .pev2Copula                 EV copula probability direct computation
#  .pevContourSlider           Interactive contour plots of EV probability
#  .pevPerspSlider             Interactive perspective plots of EV probability
# FUNCTION:                  EXTREME VALUE COPULAE DENSITY:
#  devCopula                  Computes extreme value copula density
#  devSlider                  Displays interactively plots of density
#  .dev1Copula                 EV copula density via dependence function
#  .dev2Copula                 EV copula density direct computation
#  .devContourSlider           Interactive contour plots of EV density
#  .devPerspSlider             Interactive perspective plots of EV density
################################################################################


################################################################################
# FUNCTION:                  EXTREME VALUE COPULAE RANDOM VARIATES:
#  revCopula                  Generates extreme value copula random variates
#  revSlider                  Displays interactively plots of random variates


revCopula <-
    function(n, param = NULL, type = evList())
{
    # Default Settings:
    subintervals = 100
    u = runif(n)

    # Match Arguments:
    type = match.arg(type)

    # Check Parameters:
    if (is.null(param)) param = evParam(type)$param

    # Random Variates:
    q = runif(n)
    v = u
    Y = seq(0, 1, length = subintervals)
    for (i in 1:n) {
        U = rep(u[i], times = subintervals)
        C.uv = pevCopula(u = U, v = Y, param, type) / U
        x = log(U)/log(U*Y)
        A = Afunc(x, param, type)
        Aderiv = .AfuncFirstDer(x, param, type)
        X = C.uv * (A + Aderiv * log(Y)/log(U*Y))
        v[i] = approx(X, Y, xout = q[i])$y
    }
    ans = cbind(u = u, v = v)

    # Add Control List:
    control = list(param = param, copula = "ev", type = type)
    attr(ans, "control")<-unlist(control)

    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


revSlider <- 
    function(B = 10)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Displays interactively perspective plots of random variates

    # FUNCTION:

    # Graphic Frame:
    par(mfrow = c(1, 1))

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

        # Plot:
        R = revCopula(N, param = param, type = type)
        plot(R, pch = 19, col = "steelblue")
        grid()
        title(main = Title)

        # Reset Frame:
        par(mfrow = c(1, 1))
    }

    # Open Slider Menu:
    setRmetricsOptions(.counter = 0)
    C = c("1 Gumbel: delta", "2 Galambos: delta", "3 Husler-Reis: delta",
          "4 Tawn: alpha", "... beta", "... r", "5 BB5: delta", "... theta")
    .sliderMenu(refresh.code,
        names = c("Copula", "N", C),
                         #   gumbel galamb h.r  tawn-tawn-tawn  bb5-bb5
        minima =      c(1, 100,   1,    0,   0,    0,   0,   1,   0,  1),
        maxima =      c(5,5000,   B,    B,   B,    1,   1,   B,   B,  B),
        resolutions = c(1, 100, .05,  .05, .05,  .01, .01,  .1,  .1, .1),
        starts =      c(1, 100,   2,    1,   1,   .5,  .5,   2,   1,  2))
}


################################################################################
# FUNCTION:                  EXTREME VALUE COPULAE PROBABILIY:
#  pevCopula                  Computes extreme value copula probability
#  pevSlider                  Displays interactively plots of probability
#  .pev1Copula                 EV copula probability via dependence function
#  .pev2Copula                 EV copula probability direct computation
#  .pevContourSlider           Interactive contour plots of EV probability
#  .pevPerspSlider             Interactive perspective plots of EV probability


pevCopula <-
    function(u = 0.5, v = u, param = NULL, type = evList(),
    output = c("vector", "list"), alternative = FALSE )
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Computes extreme value copula probability

    # Arguments:
    #   u, v - two numeric values or vectors of the same length at
    #       which the copula will be computed. If 'u' is a list then the
    #       the '$x' and '$y' elements will be used as 'u' and 'v'.
    #       If 'u' is a two column matrix then the first column will
    #       be used as 'u' and the the second as 'v'.
    #   param - a numeric value or vector of named parameters as
    #       required by the copula specified by the variable 'type'.
    #       If set to NULL, then the parameters will be taken as
    #       specified by the function 'evParam'.
    #   type - the type of the maximum extreme value copula. A character
    #       string selected from: "gumbel", "galambos", "husler.reiss",
    #       "tawn", or "bb5".
    #   output - a character string specifying how the output should
    #       be formatted. By default a vector of the same length as
    #       'u' and 'v'. If specified as "list" then 'u' and 'v' are
    #       expected to span a two-dimensional grid as outputted by the
    #       function 'grid2d' and the function returns a list with
    #       elements '$x', 'y', and 'z' which can be directly used
    #       for example by 2D plotting functions.
    #   alternative - Should the probability be computed alternatively
    #       in a direct way from the probability formula or by default
    #       via the dependency function?

    # Value:
    #   returns a vector or list of probabilities depending on the
    #   value of the "output" variable.

    # Example:
    #   Diagonal Value: pevCopula((0:10)/10)
    #   persp(pevCopula(u=grid2d(), output="list"), theta=-40, phi=30, xlab="x")

    # FUNCTION:

    # Select Type:
    type = match.arg(type)

    # Compute Copula:
    if (!alternative) {
        ans = .pev1Copula(u, v, param, type, output)
    } else {
        ans = .pev2Copula(u, v, param, type, output)
    }

    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


pevSlider <- 
    function(type = c("persp", "contour"), B = 10)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Displays interactively plots of probability

    # Arguments:
    #   type - a character string specifying the plot type.
    #       Either a perspective plot which is the default or
    #       a contour plot with an underlying image plot will
    #       be created.
    #   B - the maximum slider menu value when the boundary
    #       value is infinite. By default this is set to 10.

    # Match Arguments:
    type = match.arg(type)

    # Plot:
    if (type == "persp")
        .pevPerspSlider(B = B)
    if (type == "contour")
        .pevContourSlider(B = B)

    # Return Value:
    invisible()
}


# ------------------------------------------------------------------------------


.pev1Copula <- 
    function(u = 0.5, v = u, param = NULL, type = evList(),
    output = c("vector", "list") )
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Computes extreme value copula probability via dependency function

    # FUNCTION:

    # Match Arguments:
    type = match.arg(type)
    output = match.arg(output)

    # Settings:
    if (is.null(param)) {
        param = evParam(type)$param
    }
    if (is.list(u)) {
        v = u$y
        u = u$x
    }
    if (is.matrix(u)) {
        v = u[, 2]
        u = u[, 1]
    }

    # Settings:
    log.u = log(u)
    log.v = log(v)
    x = log.u/(log.u+log.v)

    # Copula Probability:
    A = Afunc(x, param = param, type = type)
    C = exp((log.u+log.v) * A)
    names(C) = NULL

    # Simulates Max function:
    C = (C + abs(C))/2

    # On Boundary:
    C[is.na(C)] = 0
    C[which(u == 0)] = 0
    C[which(u == 1)] = v[which(u == 1)]
    C[which(v == 0)] = 0
    C[which(v == 1)] = u[which(v == 1)]
    C[which(u*v == 1)] = 1
    C[which(u+v == 0)] = 0

    # Result:
    attr(C, "control") <- unlist(list(param = param, type = type))

    # As List ?
    if (output == "list") {
        N = sqrt(length(u))
        x = u[1:N]
        y = matrix(v, ncol = N)[1, ]
        C = list(x = x, y = y,  z = matrix(C, ncol = N))
    }

    # Return Value:
    C
}


# ------------------------------------------------------------------------------


.pev2Copula <- 
    function(u = 0.5, v = u, param = NULL, type = evList(),
    output = c("vector", "list") )
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Computes extreme value copula probability directly

    # FUNCTION:

    # Match Arguments:
    type = match.arg(type)
    output = match.arg(output)

    # Settings:
    if (is.null(param)) {
        param = evParam(type)$param
    }
    if (is.list(u)) {
        v = u$y
        u = u$x
    }
    if (is.matrix(u)) {
        v = u[, 2]
        u = u[, 1]
    }

    # Compute Probability:
    if (type == "gumbel") {
        alpha = param[1]
        C = exp(-((-log(u))^alpha + (-log(v))^alpha)^(1/alpha))
    }
    if (type == "galambos") {
        alpha = param[1]
        u.tilde = -log(u)
        v.tilde = -log(v)
        C = u*v*exp(((u.tilde)^(-alpha) +
            (v.tilde)^(-alpha))^(-1/alpha))
    }
    if (type == "husler.reiss") {
        alpha = param[1]
        u.tilde = -log(u)
        v.tilde = -log(v)
        C = exp(-
            u.tilde * pnorm(1/alpha + 0.5*alpha*log(u.tilde/v.tilde)) -
            v.tilde * pnorm(1/alpha + 0.5*alpha*log(v.tilde/u.tilde)) )
    }
    if (type == "tawn") {
        b = param[1]
        a = param[2]
        r = param[3]
        log.uv = log(u*v)
        t = log(u)/log.uv
        A = 1-b+(b-a)*t+(a^r*t^r+b^r*(1-t)^r)^(1/r)
        C = exp(log.uv*A)
    }
    if (type == "bb5") {
        delta = param[1]
        theta = param[2]
        u.tilde = -log(u)
        v.tilde = -log(v)
        C = exp(-(  u.tilde^theta + v.tilde^theta -
            ( u.tilde^(-theta*delta) +
              v.tilde^(-theta*delta) )^(-1/delta))^(1/theta))
    }

    # Some more, yet untested and undocumented:
    if (type == "gumbelII") {
        alpha = param[1]
        C = u*v*exp(alpha*log(u)*log(v)/(log(u)+log(v)))
    }
    if (type == "marshall.olkin") {
        a = param[1]
        b = param[2]
        C = apply(cbind(v*u^(1-a), u*v^(1-b)), 1, min)
    }
    if (type == "pi" || type == "Cperp") {
        C = u*v
    }
    if (type == "m" || type == "Cplus") {
        C = apply(cbind(u, v), 1, min)
    }

    # Simulates Max function:
    C = (C + abs(C))/2

    # On Boundary:
    C[is.na(C)] = 0
    C[which(u == 0)] = 0
    C[which(u == 1)] = v[which(u == 1)]
    C[which(v == 0)] = 0
    C[which(v == 1)] = u[which(v == 1)]
    C[which(u*v == 1)] = 1
    C[which(u+v == 0)] = 0

    # Result:
    attr(C, "control") <- unlist(list(param = param, type = type))

    # As List ?
    if (output == "list") {
        N = sqrt(length(u))
        x = u[1:N]
        y = matrix(v, ncol = N)[1, ]
        C = list(x = x, y = y,  z = matrix(C, ncol = N))
    }

    # Return Value:
    C
}


# ------------------------------------------------------------------------------


.pevContourSlider <- 
    function(B = 10)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Displays interactively contour plots of probability

    #FUNCTION:

    # Graphic Frame:
    par(mfrow = c(1, 1))

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
        nlev = .sliderMenu(no = 11)
        ncol = .sliderMenu(no = 12)

        # Title:
        type = Type[Copula]
        subTitle = paste(paste(names(param) , "="), param, collapse = " | " )
        Title = paste(" ", type, "\n", subTitle)

        # Plot:
        uv = grid2d(x = (0:N)/N)
        D = .pev1Copula(u = uv, type = type, param = param, output = "list")
        image(D, col = heat.colors(ncol) )
        contour(D, nlevels = nlev, add = TRUE)
        title(main = Title)

        # Reset Frame:
        par(mfrow = c(1, 1))
    }

    # Open Slider Menu:
    setRmetricsOptions(.counter = 0)
    C = c("1 Gumbel: delta", "2 Galambos: delta", "3 Husler-Reis: delta",
          "4 Tawn: alpha", "... beta", "... r", "5 BB5: delta", "... theta",
          "Plot - levels", "... colors")
    .sliderMenu(refresh.code,
        names = c("Copula","N", C), #gal   hr   tawn          bb5    nlev  ncol
        minima =      c(1,  10,   1,   0,   0,   0,   0,  1,  0,  1,   5,   12),
        maxima =      c(5, 100,   B,   B,   B,   1,   1,  B,  B,  B, 100,  256),
        resolutions = c(1,   1, .05, .05, .05, .01, .01, .1, .1, .1,   5,    1),
        starts =      c(1,  25,   2,   1,   1,  .5,  .5,  2,  1,  2,  10,   12))
}


# ------------------------------------------------------------------------------


.pevPerspSlider <- 
    function(B = 10)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Displays interactively perspective plots of probability

    #FUNCTION:

    # Graphic Frame:
    par(mfrow = c(1, 1))

    # Internal Function:
    refresh.code = function(...)
    {
        # Startup Counter:
        .counter <- getRmetricsOptions(".counter") + 1
        setRmetricsOptions(.counter = .counter)
        if (.counter < 12) return ()

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
        theta = .sliderMenu(no = 11)
        phi = .sliderMenu(no = 12)

        # Title:
        type = Type[Copula]
        subTitle = paste(paste(names(param) , "="), param, collapse = " | " )
        Title = paste(" ", type, "\n", subTitle)

        # Plot:
        uv = grid2d(x = (0:N)/N)
        D =  .pev1Copula(u = uv, type = type, param = param, output = "list")
        #D2 = .pev2Copula(u = uv, type = type, param = param, output = "list")
        persp(D, theta = theta, phi = phi, col = "steelblue", shade = 0.5,
            ticktype = "detailed", cex = 0.5)
        title(main = Title)

        # Reset Frame:
        par(mfrow = c(1, 1))
    }

    # Open Slider Menu:
    setRmetricsOptions(.counter = 0)
    C = c("1 Gumbel: delta", "2 Galambos: delta", "3 Husler-Reis: delta",
          "4 Tawn: alpha", "... beta", "... r", "5 BB5: delta", "... theta",
          "Plot - theta", "... phi")
    .sliderMenu(refresh.code,
        names = c("Copula", "N", C), #gal  hr  tawn          bb5    theta  phi
        minima =      c(1,  10,   1,   0,   0,   0,   0,  1,  0,  1, -180,   0),
        maxima =      c(5, 100,   B,   B,   B,   1,   1,  B,  B,  B,  180, 360),
        resolutions = c(1,   1, .05, .05, .05, .01, .01, .1, .1, .1,    1,   1),
        starts =      c(1,  25,   2,   1,   1,  .5,  .5,  2,  1,  2,  -40,  30))
}


################################################################################
# FUNCTION:                  EXTREME VALUE COPULAE DENSITY:
#  devCopula                  Computes extreme value copula density
#  devSlider                  Displays interactively plots of density
#  .dev1Copula                 EV copula density via dependence function
#  .dev2Copula                 EV copula density direct computation
#  .devContourSlider           Interactive contour plots of EV density
#  .devPerspSlider             Interactive perspective plots of EV density


devCopula <- 
    function(u = 0.5, v = u, param = NULL, type = evList(),
    output = c("vector", "list"), alternative = FALSE )
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Computes extreme value copula density from dependence function

    # Arguments:
    #   u, v - two numeric values or vectors of the same length at
    #       which the copula will be computed. If 'u' is a list then the
    #       the '$x' and '$y' elements will be used as 'u' and 'v'.
    #       If 'u' is a two column matrix then the first column will
    #       be used as 'u' and the the second as 'v'.
    #   param - a numeric value or vector of named parameters as
    #       required by the copula specified by the variable 'type'.
    #       If set to NULL, then the parameters will be taken as
    #       specified by the function 'evParam'.
    #   type - the type of the maximum extreme value copula. A character
    #       string selected from: "gumbel", "galambos", "husler.reiss",
    #       "tawn", or "bb5".
    #   output - a character string specifying how the output should
    #       be formatted. By default a vector of the same length as
    #       'u' and 'v'. If specified as "list" then 'u' and 'v' are
    #       expected to span a two-dimensional grid as outputted by the
    #       function 'grid2d' and the function returns a list with
    #       elements '$x', 'y', and 'z' which can be directly used
    #       for example by 2D plotting functions.
    #   alternative - Should the density be computed alternatively
    #       in a direct way from the probability formula or by default
    #       via the dependency function?

    # Value:
    #   returns a vector or list of density values depending on the
    #   value of the "output" variable.

    # Example:
    #   Diagonal Value: devCopula((0:10)/10)
    #   persp(devCopula(u=grid2d(), output="list"), theta=-40, phi=30, xlab="x")

    # FUNCTION:

    # Match Arguments:
    type = match.arg(type)
    output = match.arg(output)

    # Copula Density:
    if (alternative) {
        ans = .dev2Copula(u, v, param, type, output)
    } else {
        ans = .dev1Copula(u, v, param, type, output)
    }

    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


devSlider =
    function(type = c("persp", "contour"), B = 10)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Displays interactively plots of probability

    # Arguments:
    #   type - a character string specifying the plot type.
    #       Either a perspective plot which is the default or
    #       a contour plot with an underlying image plot will
    #       be created.
    #   B - the maximum slider menu value when the boundary
    #       value is infinite. By default this is set to 10.

    # Match Arguments:
    type = match.arg(type)

    # Plot:
    if (type == "persp")
        .devPerspSlider(B = B)
    if (type == "contour")
        .devContourSlider(B = B)

    # Return Value:
    invisible()
}


# ------------------------------------------------------------------------------


.dev1Copula <- 
    function(u = 0.5, v = u, param = NULL, type = evList(),
    output = c("vector", "list") )
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Computes extreme value copula density from dependence function

    # Example:
    #   Diagonal Value: devCopula((0:10)/10)
    #   persp(devCopula(u=grid2d(), output="list"), theta=-40, phi=30, xlab="x")

    # FUNCTION:

    # Match Arguments:
    type = match.arg(type)
    output = match.arg(output)

    # Settings:
    if (is.null(param)) {
        param = evParam(type)$param
    }
    if (is.list(u)) {
        v = u$y
        u = u$x
    }
    if (is.matrix(u)) {
        v = u[, 2]
        u = u[, 1]
    }

    # Settings for Maple Output:
    Pi = pi
    ln = function(x) { log(x) }
    erf = function (x) { 2*pnorm(sqrt(2)*x)-1 }

    # Further Settings:
    log.u = log(u)
    log.v = log(v)
    x = log.u/(log.u+log.v)
    y = log.v/(log.u+log.v)

    # Copula Probability:
    A = Afunc(x, param = param, type = type)
    A1 = .AfuncFirstDer(x, param = param, type = type)
    A2 = .AfuncSecondDer(x, param = param, type = type)

    # Prefactor:
    P = pevCopula(u, v, param = param, type = type) / (u*v)
    c.uv = P * (( -x*y/(log.u+log.v))*A2 + (A+y*A1)*(A-x*A1) )
    c.uv[which(u*v == 0 | u*v == 1)] = 0

    # Result:
    attr(c.uv, "control") <- unlist(list(param = param, type = type))

    # As List ?
    if (output == "list") {
        N = sqrt(length(u))
        x = u[1:N]
        y = matrix(v, ncol = N)[1, ]
        c.uv = list(x = x, y = y,  z = matrix(c.uv, ncol = N))
    }

    # Return Value:
    c.uv
}


# ------------------------------------------------------------------------------


.dev2Copula <- 
    function(u = 0.5, v = u, param = NULL, type = evList(),
    output = c("vector", "list") )
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Computes extreme value copula density directly

    # Details:
    #   List - 9 Types:
    #   pi[Cperp], gumbel, gumbelII, galambos, husler.reiss,
    #   tawn, bb5, marshall.olkin, m[Cplus]

    # References:
    #   Carmona, Evanesce

    # FUNCTION:

    # Match Arguments:
    type = match.arg(type)
    output = match.arg(output)

    # Settings:
    if (is.null(param)) {
        param = evParam(type)$param
    }
    if (is.list(u)) {
        v = u$y
        u = u$x
    }
    if (is.matrix(u)) {
        v = u[, 2]
        u = u[, 1]
    }

    # Settings:
    if (is.null(param)) param = evParam[[type]]
    Pi = pi
    ln = function(x) { log(x) }
    erf = function (x) { 2*pnorm(sqrt(2)*x)-1 }

    # Compute Probability:
    if (type == "gumbel") {
        alpha = param[1]
        # Maple Generated Output:
        c.uv =
        -((-ln(u))^alpha+(-ln(v))^alpha)^(1/alpha)*(-ln(v))^alpha/v/ln(v)/(
        (-ln(u))^alpha+(-ln(v))^alpha)^2*(-ln(u))^alpha/u/ln(u)*exp(-((-ln(
        u))^alpha+(-ln(v))^alpha)^(1/alpha))+((-ln(u))^alpha+(-ln(v))^alpha
        )^(1/alpha)*(-ln(u))^alpha/u/ln(u)/((-ln(u))^alpha+(-ln(v))^alpha)^
        2*exp(-((-ln(u))^alpha+(-ln(v))^alpha)^(1/alpha))*(-ln(v))^alpha*
        alpha/v/ln(v)+(((-ln(u))^alpha+(-ln(v))^alpha)^(1/alpha))^2*(-ln(u)
        )^alpha/u/ln(u)/((-ln(u))^alpha+(-ln(v))^alpha)^2*(-ln(v))^alpha/v/
        ln(v)*exp(-((-ln(u))^alpha+(-ln(v))^alpha)^(1/alpha))
    }
    if (type == "galambos") {
        alpha = param[1]
        # Maple Generated Output:
        c.uv =
        exp(((-ln(u))^(-alpha)+(-ln(v))^(-alpha))^(-1/alpha))+((-ln(u))^(-
        alpha)+(-ln(v))^(-alpha))^(-1/alpha)*(-ln(v))^(-alpha)/ln(v)/((-ln(
        u))^(-alpha)+(-ln(v))^(-alpha))*exp(((-ln(u))^(-alpha)+(-ln(v))^(-
        alpha))^(-1/alpha))+((-ln(u))^(-alpha)+(-ln(v))^(-alpha))^(-1/alpha
        )*(-ln(u))^(-alpha)/ln(u)/((-ln(u))^(-alpha)+(-ln(v))^(-alpha))*exp(
        ((-ln(u))^(-alpha)+(-ln(v))^(-alpha))^(-1/alpha))+((-ln(u))^(-
        alpha)+(-ln(v))^(-alpha))^(-1/alpha)*(-ln(v))^(-alpha)/ln(v)/((-ln(
        u))^(-alpha)+(-ln(v))^(-alpha))^2*(-ln(u))^(-alpha)/ln(u)*exp(((-ln
        (u))^(-alpha)+(-ln(v))^(-alpha))^(-1/alpha))+((-ln(u))^(-alpha)+(-
        ln(v))^(-alpha))^(-1/alpha)*(-ln(u))^(-alpha)/ln(u)/((-ln(u))^(-
        alpha)+(-ln(v))^(-alpha))^2*exp(((-ln(u))^(-alpha)+(-ln(v))^(-alpha
        ))^(-1/alpha))*(-ln(v))^(-alpha)*alpha/ln(v)+(((-ln(u))^(-alpha)+(-
        ln(v))^(-alpha))^(-1/alpha))^2*(-ln(u))^(-alpha)/ln(u)/((-ln(u))^(-
        alpha)+(-ln(v))^(-alpha))^2*(-ln(v))^(-alpha)/ln(v)*exp(((-ln(u))^(
        -alpha)+(-ln(v))^(-alpha))^(-1/alpha))
    }
    if (type == "husler.reiss") {
        # Maple Generated Output:
        c.uv =
        (-.2500000000/u/Pi^(1/2)*exp(-1/2*(1/alpha+.5*alpha*ln(ln(u)/ln(v))
        )^2)*alpha/v/ln(v)*2^(1/2)+.1250000000/Pi^(1/2)*(1/alpha+.5*alpha*
        ln(ln(u)/ln(v)))*alpha^2/v/ln(v)*exp(-1/2*(1/alpha+.5*alpha*ln(ln(u
        )/ln(v)))^2)/u*2^(1/2)-.2500000000/v/Pi^(1/2)*exp(-1/2*(1/alpha+.5*
        alpha*ln(ln(v)/ln(u)))^2)*alpha/u/ln(u)*2^(1/2)+.1250000000/Pi^(1/2
        )*(1/alpha+.5*alpha*ln(ln(v)/ln(u)))*alpha^2/v*exp(-1/2*(1/alpha+.5
        *alpha*ln(ln(v)/ln(u)))^2)/u/ln(u)*2^(1/2))*exp(.5*ln(u)*(erf(1/2*(
        1/alpha+.5*alpha*ln(ln(u)/ln(v)))*2^(1/2))+1)+.5*ln(v)*(erf(1/2*(1/
        alpha+.5*alpha*ln(ln(v)/ln(u)))*2^(1/2))+1))+(.5/u*(erf(1/2*(1/
        alpha+.5*alpha*ln(ln(u)/ln(v)))*2^(1/2))+1)+.2500000000/Pi^(1/2)*
        exp(-1/2*(1/alpha+.5*alpha*ln(ln(u)/ln(v)))^2)*alpha/u*2^(1/2)-.25*
        ln(v)/Pi^(1/2)*exp(-1/2*(1/alpha+.5*alpha*ln(ln(v)/ln(u)))^2)*
        alpha/u/ln(u)*2^(1/2))*(-.2500000000*ln(u)/Pi^(1/2)*exp(-1/2*(1/
        alpha+.5*alpha*ln(ln(u)/ln(v)))^2)*alpha/v/ln(v)*2^(1/2)+.5/v*(erf(
        1/2*(1/alpha+.5*alpha*ln(ln(v)/ln(u)))*2^(1/2))+1)+.2500000000/Pi^(
        1/2)*exp(-1/2*(1/alpha+.5*alpha*ln(ln(v)/ln(u)))^2)*alpha/v*2^(1/2)
        )*exp(.5*ln(u)*(erf(1/2*(1/alpha+.5*alpha*ln(ln(u)/ln(v)))*2^(1/2))
        +1)+.5*ln(v)*(erf(1/2*(1/alpha+.5*alpha*ln(ln(v)/ln(u)))*2^(1/2))+1
        ))
    }
    if (type == "tawn") {
        # 0 <= alpha, beta <= 1, 1 <= r < Inf
        b = param[1]
        a = param[2]
        r = param[3]
        # Maple Generated Output:
        c.uv =
        (-(b-a)/u/ln(u*v)^2/v+2*(b-a)*ln(u)/ln(u*v)^3/u/v+(a^r*(ln(u)/ln(u*
        v))^r+b^r*(1-ln(u)/ln(u*v))^r)^(1/r)/r^2*(-a^r*(ln(u)/ln(u*v))^r*r/
        ln(u*v)/v+b^r*(1-ln(u)/ln(u*v))^r*r*ln(u)/ln(u*v)^2/v/(1-ln(u)/ln(u
        *v)))/(a^r*(ln(u)/ln(u*v))^r+b^r*(1-ln(u)/ln(u*v))^r)^2*(a^r*(ln(u)
        /ln(u*v))^r*r*(1/u/ln(u*v)-ln(u)/ln(u*v)^2/u)/ln(u)*ln(u*v)+b^r*(1-
        ln(u)/ln(u*v))^r*r*(-1/u/ln(u*v)+ln(u)/ln(u*v)^2/u)/(1-ln(u)/ln(u*v
        )))+(a^r*(ln(u)/ln(u*v))^r+b^r*(1-ln(u)/ln(u*v))^r)^(1/r)/r*(-a^r*(
        ln(u)/ln(u*v))^r*r^2/v*(1/u/ln(u*v)-ln(u)/ln(u*v)^2/u)/ln(u)+a^r*(
        ln(u)/ln(u*v))^r*r*(-1/u/ln(u*v)^2/v+2*ln(u)/ln(u*v)^3/u/v)/ln(u)*
        ln(u*v)+a^r*(ln(u)/ln(u*v))^r*r*(1/u/ln(u*v)-ln(u)/ln(u*v)^2/u)/ln(
        u)/v+b^r*(1-ln(u)/ln(u*v))^r*r^2*ln(u)/ln(u*v)^2/v/(1-ln(u)/ln(u*v)
        )^2*(-1/u/ln(u*v)+ln(u)/ln(u*v)^2/u)+b^r*(1-ln(u)/ln(u*v))^r*r*(1/u
        /ln(u*v)^2/v-2*ln(u)/ln(u*v)^3/u/v)/(1-ln(u)/ln(u*v))-b^r*(1-ln(u)/
        ln(u*v))^r*r*(-1/u/ln(u*v)+ln(u)/ln(u*v)^2/u)/(1-ln(u)/ln(u*v))^2*
        ln(u)/ln(u*v)^2/v)/(a^r*(ln(u)/ln(u*v))^r+b^r*(1-ln(u)/ln(u*v))^r)-
        (a^r*(ln(u)/ln(u*v))^r+b^r*(1-ln(u)/ln(u*v))^r)^(1/r)/r*(a^r*(ln(u)
        /ln(u*v))^r*r*(1/u/ln(u*v)-ln(u)/ln(u*v)^2/u)/ln(u)*ln(u*v)+b^r*(1-
        ln(u)/ln(u*v))^r*r*(-1/u/ln(u*v)+ln(u)/ln(u*v)^2/u)/(1-ln(u)/ln(u*v
        )))/(a^r*(ln(u)/ln(u*v))^r+b^r*(1-ln(u)/ln(u*v))^r)^2*(-a^r*(ln(u)/
        ln(u*v))^r*r/ln(u*v)/v+b^r*(1-ln(u)/ln(u*v))^r*r*ln(u)/ln(u*v)^2/v/
        (1-ln(u)/ln(u*v))))*exp(ln(u*v)-b+(b-a)*ln(u)/ln(u*v)+(a^r*(ln(u)/
        ln(u*v))^r+b^r*(1-ln(u)/ln(u*v))^r)^(1/r))+(1/u+(b-a)/u/ln(u*v)-(b-
        a)*ln(u)/ln(u*v)^2/u+(a^r*(ln(u)/ln(u*v))^r+b^r*(1-ln(u)/ln(u*v))^r
        )^(1/r)/r*(a^r*(ln(u)/ln(u*v))^r*r*(1/u/ln(u*v)-ln(u)/ln(u*v)^2/u)/
        ln(u)*ln(u*v)+b^r*(1-ln(u)/ln(u*v))^r*r*(-1/u/ln(u*v)+ln(u)/ln(u*v)
        ^2/u)/(1-ln(u)/ln(u*v)))/(a^r*(ln(u)/ln(u*v))^r+b^r*(1-ln(u)/ln(u*v
        ))^r))*(1/v-(b-a)*ln(u)/ln(u*v)^2/v+(a^r*(ln(u)/ln(u*v))^r+b^r*(1-
        ln(u)/ln(u*v))^r)^(1/r)/r*(-a^r*(ln(u)/ln(u*v))^r*r/ln(u*v)/v+b^r*(
        1-ln(u)/ln(u*v))^r*r*ln(u)/ln(u*v)^2/v/(1-ln(u)/ln(u*v)))/(a^r*(ln(
        u)/ln(u*v))^r+b^r*(1-ln(u)/ln(u*v))^r))*exp(ln(u*v)-b+(b-a)*ln(u)/
        ln(u*v)+(a^r*(ln(u)/ln(u*v))^r+b^r*(1-ln(u)/ln(u*v))^r)^(1/r))
    }
    if (type == "bb5") {
        # delta > 0, theta >= 1
        delta = param[1]
        theta = param[2]
        # Maple Generated Output:
        c.uv =
        -((-ln(u))^theta+(-ln(v))^theta-((-ln(u))^(-theta*delta)+(-ln(v))^(
        -theta*delta))^(-1/delta))^(1/theta)/theta^2*((-ln(v))^theta*theta/
        v/ln(v)-((-ln(u))^(-theta*delta)+(-ln(v))^(-theta*delta))^(-1/delta
        )*(-ln(v))^(-theta*delta)*theta/v/ln(v)/((-ln(u))^(-theta*delta)+(-
        ln(v))^(-theta*delta)))/((-ln(u))^theta+(-ln(v))^theta-((-ln(u))^(-
        theta*delta)+(-ln(v))^(-theta*delta))^(-1/delta))^2*((-ln(u))^theta
        *theta/u/ln(u)-((-ln(u))^(-theta*delta)+(-ln(v))^(-theta*delta))^(-
        1/delta)*(-ln(u))^(-theta*delta)*theta/u/ln(u)/((-ln(u))^(-theta*
        delta)+(-ln(v))^(-theta*delta)))*exp(-((-ln(u))^theta+(-ln(v))^
        theta-((-ln(u))^(-theta*delta)+(-ln(v))^(-theta*delta))^(-1/delta))
        ^(1/theta))-((-ln(u))^theta+(-ln(v))^theta-((-ln(u))^(-theta*delta)
        +(-ln(v))^(-theta*delta))^(-1/delta))^(1/theta)/theta*(-((-ln(u))^(
        -theta*delta)+(-ln(v))^(-theta*delta))^(-1/delta)*(-ln(v))^(-theta*
        delta)*theta^2/v/ln(v)/((-ln(u))^(-theta*delta)+(-ln(v))^(-theta*
        delta))^2*(-ln(u))^(-theta*delta)/u/ln(u)-((-ln(u))^(-theta*delta)+
        (-ln(v))^(-theta*delta))^(-1/delta)*(-ln(u))^(-theta*delta)*theta^2
        /u/ln(u)/((-ln(u))^(-theta*delta)+(-ln(v))^(-theta*delta))^2*(-ln(v
        ))^(-theta*delta)*delta/v/ln(v))/((-ln(u))^theta+(-ln(v))^theta-((-
        ln(u))^(-theta*delta)+(-ln(v))^(-theta*delta))^(-1/delta))*exp(-((-
        ln(u))^theta+(-ln(v))^theta-((-ln(u))^(-theta*delta)+(-ln(v))^(-
        theta*delta))^(-1/delta))^(1/theta))+((-ln(u))^theta+(-ln(v))^theta
        -((-ln(u))^(-theta*delta)+(-ln(v))^(-theta*delta))^(-1/delta))^(1/
        theta)/theta*((-ln(u))^theta*theta/u/ln(u)-((-ln(u))^(-theta*delta)
        +(-ln(v))^(-theta*delta))^(-1/delta)*(-ln(u))^(-theta*delta)*theta/
        u/ln(u)/((-ln(u))^(-theta*delta)+(-ln(v))^(-theta*delta)))/((-ln(u)
        )^theta+(-ln(v))^theta-((-ln(u))^(-theta*delta)+(-ln(v))^(-theta*
        delta))^(-1/delta))^2*exp(-((-ln(u))^theta+(-ln(v))^theta-((-ln(u))
        ^(-theta*delta)+(-ln(v))^(-theta*delta))^(-1/delta))^(1/theta))*((-
        ln(v))^theta*theta/v/ln(v)-((-ln(u))^(-theta*delta)+(-ln(v))^(-
        theta*delta))^(-1/delta)*(-ln(v))^(-theta*delta)*theta/v/ln(v)/((-
        ln(u))^(-theta*delta)+(-ln(v))^(-theta*delta)))+(((-ln(u))^theta+(-
        ln(v))^theta-((-ln(u))^(-theta*delta)+(-ln(v))^(-theta*delta))^(-1/
        delta))^(1/theta))^2/theta^2*((-ln(u))^theta*theta/u/ln(u)-((-ln(u)
        )^(-theta*delta)+(-ln(v))^(-theta*delta))^(-1/delta)*(-ln(u))^(-
        theta*delta)*theta/u/ln(u)/((-ln(u))^(-theta*delta)+(-ln(v))^(-
        theta*delta)))/((-ln(u))^theta+(-ln(v))^theta-((-ln(u))^(-theta*
        delta)+(-ln(v))^(-theta*delta))^(-1/delta))^2*((-ln(v))^theta*theta
        /v/ln(v)-((-ln(u))^(-theta*delta)+(-ln(v))^(-theta*delta))^(-1/
        delta)*(-ln(v))^(-theta*delta)*theta/v/ln(v)/((-ln(u))^(-theta*
        delta)+(-ln(v))^(-theta*delta)))*exp(-((-ln(u))^theta+(-ln(v))^
        theta-((-ln(u))^(-theta*delta)+(-ln(v))^(-theta*delta))^(-1/delta))
        ^(1/theta))
    }

    # Result:
    attr(c.uv, "control") <- unlist(list(param = param, type = type))

    # As List ?
    if (output[1] == "list") {
        N = sqrt(length(u))
        x = u[1:N]
        y = matrix(v, ncol = N)[1, ]
        c.uv = list(x = x, y = y,  z = matrix(c.uv, ncol = N))
    }

    # Return Value:
    c.uv
}


# ------------------------------------------------------------------------------


.devContourSlider <- 
    function(B = 10)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Displays interactively contour plots of density

    # FUNCTION:

    # Internal Function:
    refresh.code = function(...)
    {
        # Startup Counter:
        .counter <- getRmetricsOptions(".counter") + 1
        setRmetricsOptions(.counter = .counter)
        if (.counter < 12) return ()

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
        nlev = .sliderMenu(no = 11)
        ncol = .sliderMenu(no = 12)

        # Title:
        type = Type[Copula]
        subTitle = paste(paste(names(param) , "="), param, collapse = " | " )
        Title = paste(" ", type, "\n", subTitle)

        # Plot:
        n = N/2
        F = (2*1.0e-2)^(1/n)
        x = 0.5*F^(1:n)
        x = c(rev(x), 0.5, 1-x)
        uv = grid2d(x = (1:(N-1))/N)
        D = .dev1Copula(u = uv, type = type, param = param, output = "list")
        image(D, col = heat.colors(ncol) )
        contour(D, nlevels = nlev, add = TRUE)
        title(main = Title)

        # Reset Frame:
        par(mfrow = c(1, 1))
    }

    # Open Slider Menu:
    setRmetricsOptions(.counter = 0)
    C = c("1 Gumbel: delta", "2 Galambos: delta", "3 Husler-Reis: delta",
          "4 Tawn: alpha", "... beta", "... r", "5 BB5: delta", "... theta",
          "Plot - levels", "... colors")
    .sliderMenu(refresh.code,
        names = c("Copula","N", C), #gal   hr   tawn          bb5    nlev  ncol
        minima =      c(1,  10,   1,   0,   0,   0,   0,  1,  0,  1,   5,   12),
        maxima =      c(5, 100,   B,   B,   B,   1,   1,  B,  B,  B, 100,  256),
        resolutions = c(1,   5, .05, .05, .05, .01, .01, .1, .1, .1,   5,    1),
        starts =      c(1,  25,   2,   1,   1,  .5,  .5,  2,  1,  2,  10,   12))
}


# ------------------------------------------------------------------------------


.devPerspSlider <-
    function(B = 10)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Displays interactively contour plots of density

    #FUNCTION:

    # Internal Function:
    refresh.code = function(...)
    {
        # Startup Counter:
        .counter <- getRmetricsOptions(".counter") + 1
        setRmetricsOptions(.counter = .counter)
        if (.counter < 12) return ()

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
        theta = .sliderMenu(no = 11)
        phi = .sliderMenu(no = 12)

        # Title:
        type = Type[Copula]
        subTitle = paste(paste(names(param) , "="), param, collapse = " | " )
        Title = paste(" ", type, "\n", subTitle)

        # Plot:
        n = N/2
        F = (2*1.0e-2)^(1/n)
        x = 0.5*F^(1:n)
        x = c(rev(x), 0.5, 1-x)
        uv = grid2d(x = x)
        D = .dev1Copula(u = uv, type = type, param = param, output = "list")
        persp(D, theta = theta, phi = phi, col = "steelblue", shade = 0.5,
            ticktype = "detailed", cex = 0.5)
        title(main = Title)

        # Reset Frame:
        par(mfrow = c(1, 1))
    }

    # Open Slider Menu:
    setRmetricsOptions(.counter = 12)
    C = c("1 Gumbel: delta", "2 Galambos: delta", "3 Husler-Reis: delta",
          "4 Tawn: alpha", "... beta", "... r", "5 BB5: delta", "... theta",
          "Plot - theta", "... phi")
    .sliderMenu(refresh.code,
        names = c("Copula", "N", C), #gal  hr  tawn          bb5    theta  phi
        minima =      c(1,  10,   1,   0,   0,   0,   0,  1,  0,  1, -180,   0),
        maxima =      c(5, 100,   B,   B,   B,   1,   1,  B,  B,  B,  180, 360),
        resolutions = c(1,   5, .05, .05, .05, .01, .01, .1, .1, .1,    1,   1),
        starts =      c(1,  25,   2,   1,   1,  .5,  .5,  2,  1,  2,  -40,  30))
}


################################################################################

