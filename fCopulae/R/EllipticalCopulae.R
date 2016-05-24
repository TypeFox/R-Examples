
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
# FUNCTION:                  ELLIPTICAL COPULAE RANDOM DEVIATES:
#  rellipticalCopula          Generates elliptical copula variates
#  rellipticalSlider          Generates interactive plots of random variates
#  .rnormCopula               Generates normal copula random variate
#  .rcauchyCopula             Generates Cauchy copula random variate
#  .rtCopula                  Generates Student-t copula random variate
# FUNCTION:                  ELLIPTICAL COPULAE PROBABILITY:
#  pellipticalCopula          Computes elliptical copula probability
#  pellipticalSlider          Generates interactive plots of probability
#  .pnormCopula               Computes normal copula probability
#  .pcauchyCopula             Computes Cauchy copula probability
#  .ptCopula                  Computes Student-t copula probability
#  .pellipticalCopulaGrid     Fast equidistant grid version
#  .pellipticalCopulaDiag     Fast diagonal cross section version
#  .pellipticalPerspSlider    Interactive perspective plots of probability
#  .pellipticalContourSlider  Interactive contour plots of probability
# FUNCTION:                  ELLIPTICAL COPULAE DENSITY:
#  dellipticalCopula          Computes elliptical copula density
#  dellipticalSlider          Generates interactive plots of density
#  .dnormCopula               Computes normal copula density
#  .dcauchyCopula             Computes Cauchy copula density
#  .dtCopula                  Computes Student-t copula density
#  .dellipticalCopulaGrid     Fast grid version for elliptical copula density
#  .dellipticalPerspSlider    Interactive perspective plots of density
#  .dellipticalContourSlider  Interactive contour plots of density
################################################################################


################################################################################
# FUNCTION:                  ELLIPTICAL COPULAE RANDOM DEVIATES:
#  rellipticalCopula          Generates elliptical copula variates
#  rellipticalSlider          Generates interactive plots of random variates
#  .rnormCopula               Generates normal copula random variate
#  .rcauchyCopula             Generates Cauchy copula random variate
#  .rtCopula                  Generates Student-t copula random variate


rellipticalCopula <-
    function(n, rho = 0.75, param = NULL, type = c("norm", "cauchy", "t"))
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Computes extreme value copula probability

    # Arguments:
    #   n - number of deviates to be generated.
    #   rho - a numeric value setting the coorelation strength, ranging
    #       between minus one and one.
    #   nu - the number of degrees of freedom, only required for
    #       Student-t copulae.
    #   type - the type of the elliptical copula. Either "norm" or
    #       "t" denoting the normal or Student-t copula, respectively.
    #   output - a character string specifying how the output should
    #       be formatted. By default a vector of the same length as
    #       'u' and 'v'. If specified as "list" then 'u' and 'v' are
    #       expected to span a two-dimensional grid as outputted by the
    #       function 'grid2d' and the function returns a list with
    #       elements '$x', 'y', and 'z' which can be directly used
    #       for example by 2D plotting functions.

    # Value:
    #   returns a vector or list of probabilities depending on the
    #   value of the "output" variable.

    # Example:
    #   Diagonal Value: pnormCopula((0:10)/10)
    #   persp(pnormCopula(u = grid2d(), output = "list"))

    # FUNCTION:

    # Settings:
    type = match.arg(type)

    # Parameters:
    if (type == "t") {
        if(is.null(param)) {
            param = c(nu = 4)
        } else {
            param = c(nu = param)
        }
        names(param) = "nu"
    }

    # Copula:
    if (type == "norm")
        ans = .rnormCopula(n = n, rho = rho)
    if (type == "cauchy")
        ans = .rcauchyCopula(n = n, rho = rho)
    if (type == "t")
        ans = .rtCopula(n = n, rho = rho, nu = param)

    # Add Control Attribute:
    control = list(rho = rho, param = param, type = type)
    attr(ans, "control")<-unlist(control)

    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


rellipticalSlider <- 
    function(B = 100)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Displays interactively perspective plots of random variates

    #FUNCTION:

    # Graphic Frame:
    par(mfrow = c(1, 1))

    # Internal Function:
    refresh.code <-  function(...)
    {
        # Startup Counter:
        .counter <- getRmetricsOptions(".counter") + 1
        setRmetricsOptions(.counter = .counter)
        if (.counter < 7) return ()

        # Sliders:
        Copula = .sliderMenu(no = 1)
        N = .sliderMenu(no = 2)
        rho = .sliderMenu(no = 3)
        nu = .sliderMenu(no = 4)
        seed = .sliderMenu(no = 5)
        size = .sliderMenu(no = 6)
        col = .sliderMenu(no = 7)
        Names = c("- Normal", "- Cauchy", "- Student t")
        Type = c("norm", "cauchy", "t")
        eps = 1.0e-6
        if (rho == +1) rho = rho - eps
        if (rho == -1) rho = rho + eps

        # Tau and Rho:
        Tau = ellipticalTau(rho)
        Rho = ellipticalRho(rho)

        # Plot:
        Title = paste("Elliptical Copula No:", as.character(Copula),
            Names[Copula], "\nrho =", as.character(rho), "|")
        if (Copula == 2) Title = paste(Title, "nu =", as.character(nu), "|")
        Title = paste(Title,
            "Kendall = ", as.character(round(Tau, digits = 3)), "|",
            "Spearman = ", as.character(round(Rho, digits = 3)) )
        set.seed(seed)
        R = rellipticalCopula(n = N, rho = rho, param = nu, type = Type[Copula])
        plot(x = R[, 1], y = R[, 2], xlim = c(0, 1), ylim = c(0, 1),
            xlab = "u", ylab = "v", pch = 19, col = col, cex = size)
        title(main = Title)

        # Reset Frame:
        par(mfrow = c(1, 1))
    }

    # Open Slider Menu:
    setRmetricsOptions(.counter = 0)
    plot.names = c("Plot - size", "... color")
    .sliderMenu(refresh.code,
        names       = c("Copula",   "N", "rho", "t: nu", "seed", plot.names),
        minima      = c(       1,  1000,    -1,       1,   1000,     0,   1),
        maxima      = c(       3, 10000,    +1,       B,   9999,     1,  16),
        resolutions = c(       1,   500,  0.01,       1,      1,   0.1,   1),
        starts      = c(       1,  1000,     0,       4,   4711,   0.5,   1))
}


# ------------------------------------------------------------------------------


.rnormCopula <- 
    function(n, rho = 0.75)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Generates normal copula random variate

    # Example:
    #   UV = rnormCopula(n = 10000); plot(UV[,1], UV[,2], cex = 0.25)

    # FUNCTION:

    # Use: X = .rnorm2d(n, rho) or alternatively:
    X = fMultivar:::.rnorm2d(n = n, rho = rho)

    # Generate
    Z <- NULL
    for(i in (1:n)) Z <- rbind(Z, pnorm(X [i,]))

    # Return Value:
    Z
}


# ------------------------------------------------------------------------------


.rcauchyCopula <-
    function(n, rho = 0.75)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Generates Student-t copula random variate

    # Example:
    #   UV = rtCopula(n = 10000); plot(UV[,1], UV[,2], cex = 0.25)

    # FUNCTION:

    # Cauchy Deviates:
    Z = .rtCopula(n = n, rho = rho, nu = 1)

    # Return Value:
    Z
}


# ------------------------------------------------------------------------------


.rtCopula <- 
    function(n, rho = 0.75, nu = 4)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Generates Student-t copula random variate

    # Example:
    #   UV = rtCopula(n = 10000); plot(UV[,1], UV[,2], cex = 0.25)

    # FUNCTION:

    # Use: X = .rnorm2d(n, rho) or alternatively:
    X = rt2d(n = n, rho = rho, nu = nu)

    # Generate
    Z = NULL
    for (i in (1:n)) Z = rbind(Z, pt(X [i,], df = nu))

    # Return Value:
    Z
}


################################################################################
# FUNCTION:                  ELLIPTICAL COPULAE PROBABILITY:
#  pellipticalCopula          Computes elliptical copula probability
#  pellipticalSlider          Generates interactive plots of probability
#  .pnormCopula               Computes normal copula probability
#  .pcauchyCopula             Computes Cauchy copula probability
#  .ptCopula                  Computes Student-t copula probability
#  .pellipticalCopulaGrid     Fast equidistant grid version
#  .pellipticalCopulaDiag     Fast diagonal cross section version
#  .pellipticalPerspSlider    Interactive perspective plots of probability
#  .pellipticalContourSlider  Interactive contour plots of probability


pellipticalCopula <- 
    function(u = 0.5, v = u, rho = 0.75, param = NULL, type = ellipticalList(),
    output = c("vector", "list"), border = TRUE)
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
    #   rho - a numeric value setting the coorelation strength, ranging
    #       between minus one and one.
    #   param - distributional parameters, the number of degrees of
    #       freedom for the Student-t copulae.
    #   type - the type of the elliptical copula. Either "norm" or
    #       "t" denoting the normal or Student-t copula, respectively.
    #   output - a character string specifying how the output should
    #       be formatted. By default a vector of the same length as
    #       'u' and 'v'. If specified as "list" then 'u' and 'v' are
    #       expected to span a two-dimensional grid as outputted by the
    #       function 'grid2d' and the function returns a list with
    #       elements '$x', 'y', and 'z' which can be directly used
    #       for example by 2D plotting functions.

    # Value:
    #   returns a vector or list of probabilities depending on the
    #   value of the "output" variable.

    # FUNCTION:

    # Match Arguments:
    type <- match.arg(type)
    output <- match.arg(output)

    # Settings:
    subdivisions = 100
    if (is.list(u)) {
        v = u[[2]]
        u = u[[1]]
    }
    if (is.matrix(u)) {
        v = u[, 2]
        u = u[, 1]
    }
    if (length(u) == 1 & u[1] > 1) {
        return(.pellipticalCopulaGrid(N = u, rho, param, type, border = border))
    }

    # Parameters:
    if (type == "t") if (is.null(param)) param = c(nu = 4)
    if (type == "kotz") if (is.null(param)) param = c(r = 1)
    if (type == "epower") if (is.null(param)) param = c(r = 1, s = 1)

    # Specical Copulae:
    if (type == "norm") {
        if (rho == -1) {
            ans <- pfrechetCopula(u = u, v = v, type = "m", output = output)
            return(ans)
        } else if (rho == +1) {
            ans <- pfrechetCopula(u = u, v = v, type = "w", output = output)
            return(ans)
        } else {
            ans = .pnormCopula(u = u, v = v, rho = rho, output = output)
            return(ans)
        }
    } else if (type == "cauchy") {
        ans <- .pcauchyCopula(u = u, v = v, rho = rho, output = output)
        return(ans)
    } else if (type == "t") {
        if (is.null(param)) param = 4
        ans <- .ptCopula(u = u, v = v, rho = rho, nu = param, output = output)
        return(ans)
    }

    # The remaining Copulae - Compute Density on Regular Grid:
    N = subdivisions
    x = (0:N)/N
    c.uv = .dellipticalCopulaGrid(N = N, rho, param, type, border = TRUE)
    c.uv$z[is.na(c.uv$z)] = 0

    # Integrate to get Probability:
    C.uv = 0*c.uv$z
    for (i in 1:(N+1)) {
        D = matrix(rep(0, times = (N+1)^2), ncol = N+1)
        for (j in 1:i) {
            D[1:i, j] = 1
            C.uv[i,j] = C.uv[j,i] = sum(D*c.uv$z)
        }
    }
    C.uv = C.uv/N^2

    # Take care about the Boundary on the Unit Square:
    C.uv[1, ] = C.uv[, 1] = 0
    C.uv[N+1, ] = C.uv[, N+1] = c.uv$x

    # Interpolate for the desired Values on the grid:
    U0 = trunc(u*N)
    V0 = trunc(v*N)
    P = (u - U0/N)
    Q = (v - V0/N)
    U0 = U0 + 1
    U1 = U0 + 1
    V0 = V0 + 1
    V1 = V0 + 1
    C.vec = rep(NA, times = length(u))
    for ( i in 1:length(u) ) {
        p = P[i]
        q = Q[i]
        if (p == 0 & q == 0) {
            C.vec[i] = C.uv[U0[i], V0[i]]
        } else if (p == 0 & q > 0) {
            C.vec[i] = (1-q)*C.uv[U0[i], V0[i]] + q*C.uv[U0[i], V1[i]]
        } else if (p > 0 & q == 0) {
            C.vec[i] = (1-p)*C.uv[U0[i], V0[i]] + p*C.uv[U1[i], V0[i]]
        } else {
            C.vec[i] = (1-p)*(1-q)*C.uv[U0[i], V0[i]] +
                p*(1-q)*C.uv[U1[i], V0[i]] + (1-p)*q*C.uv[U0[i], V1[i]] +
                p*q*C.uv[U1[i], V1[i]]
        }
    }
    C.uv = round(C.vec, digits = 3)
    attr(C.uv, "control") <- c(rho = rho)

    # As List ?
    if (output == "list") {
        N = sqrt(length(u))
        x = u[1:N]
        y = matrix(v, ncol = N)[1, ]
        names(x) = names(y) = NULL
        C.uv = list(x = x, y = y,  z = matrix(C.uv, ncol = N))
    }

    # Return Value:
    C.uv
}


# ------------------------------------------------------------------------------


pellipticalSlider <-
    function(type = c("persp", "contour"), B = 20)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Displays interactively plots of probability

    # Description:
    #   Displays interactively plots of probability

    # Arguments:
    #   type - a character string specifying the plot type.
    #       Either a perspective plot which is the default or
    #       a contour plot with an underlying image plot will
    #       be created.
    #   B - the maximum slider menu value when the boundary
    #       value is infinite. By default this is set to 10.

    # FUNCTION:

    # Settings:
    type = match.arg(type)

    # Plot:
    if (type == "persp")
        .pellipticalPerspSlider(B = B)
    if (type == "contour")
        .pellipticalContourSlider(B = B)

    # Return Value:
    invisible()
}


# ------------------------------------------------------------------------------


.pnormCopula <-
    function(u = 0.5, v = u, rho = 0.75, output = c("vector", "list") )
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Computes normal copula probability

    # Arguments:
    #   see function 'pellipticalCopula'

    # FUNCTION:

    # Type:
    output = match.arg(output)

    # Settings:
    type = "norm"
    if (is.list(u)) {
        v = u[[2]]
        u = u[[1]]
    }
    if (is.matrix(u)) {
        v = u[, 2]
        u = u[, 1]
    }

    # Copula Probability:
    C.uv = pnorm2d(qnorm(u), qnorm(v), rho = rho)
    names(C.uv) = NULL

    # Simulates Max function:
    C.uv = (C.uv + abs(C.uv))/2

    # On Boundary:
    C.uv[is.na(C.uv)] = 0
    C.uv[which(u == 0)] = 0
    C.uv[which(u == 1)] = v[which(u == 1)]
    C.uv[which(v == 0)] = 0
    C.uv[which(v == 1)] = u[which(v == 1)]
    C.uv[which(u*v == 1)] = 1
    C.uv[which(u+v == 0)] = 0

    # Result:
    attr(C.uv, "control") <- c(rho = rho)

    # As List ?
    if (output == "list") {
        N = sqrt(length(u))
        x = u[1:N]
        y = matrix(v, ncol = N)[1, ]
        C.uv = list(x = x, y = y,  z = matrix(C.uv, ncol = N))
    }

    # Return Value:
    C.uv
}


# ------------------------------------------------------------------------------


.pcauchyCopula <- 
    function(u = 0.5, v = u, rho = 0.75, output = c("vector", "list") )
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Computes Student-t copula probability

    # Arguments:
    #   see function 'pellipticalCopula'

    # FUNCTION:

    # Cauchy Probability:
    C.uv <- .ptCopula(u = u, v = v, rho = rho, nu = 1, output = output)
    attr(C.uv, "control") <- c(rho = rho)

    # Return Value:
    C.uv
}


# ------------------------------------------------------------------------------


.ptCopula <-
    function(u = 0.5, v = u, rho = 0.75, nu = 4, output = c("vector", "list") )
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Computes Student-t copula probability

    # Arguments:
    #   see function 'pellipticalCopula'

    # FUNCTION:

    # Match Arguments:
    output <- match.arg(output)

    # Settings:
    type = "t"
    if (is.list(u)) {
        v = u[[2]]
        u = u[[1]]
    }
    if (is.matrix(u)) {
        v = u[, 2]
        u = u[, 1]
    }

    # Copula Probability:
    C.uv <- pt2d(qt(u, df = nu), qt(v, df = nu), rho = rho, nu = nu)
    names(C.uv) = NULL

    # Simulates Max function:
    C.uv = (C.uv + abs(C.uv))/2

    # On Boundary:
    C.uv[is.na(C.uv)] = 0
    C.uv[which(u == 0)] = 0
    C.uv[which(u == 1)] = v[which(u == 1)]
    C.uv[which(v == 0)] = 0
    C.uv[which(v == 1)] = u[which(v == 1)]
    C.uv[which(u*v == 1)] = 1
    C.uv[which(u+v == 0)] = 0

    # Result:
    attr(C.uv, "control") <- c(rho = rho, nu = nu)

    # As List ?
    if (output == "list") {
        N = sqrt(length(u))
        x = u[1:N]
        y = matrix(v, ncol = N)[1, ]
        C.uv = list(x = x, y = y,  z = matrix(C.uv, ncol = N))
    }

    # Return Value:
    C.uv
}


# ------------------------------------------------------------------------------


.pellipticalCopulaGrid <-
    function(N, rho = 0.75, param = NULL, type = ellipticalList(), border = TRUE)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Computes elliptical copula probability on a 2d grid

    # Arguments:
    #   see function pellipticalCopula()

    # FUNCTION:

    # Settings:
    U = (0:N)/N
    V = (1:(N-1))/N

    # Compute Density on Regular Grid:
    c.uv = .dellipticalCopulaGrid(N, rho, param, type, border = TRUE)
    c.uv$z[is.na(c.uv$z)] = 0

    # Integrate to get Probability:
    if (TRUE) {
        C.uv = 0*c.uv$z
        for (i in 1:(N+1)) {
            for (j in 1:i) {
                C.uv[i,j] = C.uv[j,i] = sum(c.uv$z[1:i, 1:j])
            }
        }
        C.uv = C.uv/N^2
    }

    if (FALSE) {
        # This is much slower !
        IJ = grid2d(1:(N+1))
        X = cbind(IJ$x, IJ$y)
        fun =  function(X, C) sum(C[1:X[1], 1:X[2]])
        C.uv = apply(X, MARGIN=1, FUN = fun, C = c.uv$z)
        C.uv = matrix(C.uv, byrow = TRUE, ncol = N+1) / N^2
    }

    # Probability - Take care about the Boundary on the Unit Square:
    C.uv[1, ] = C.uv[, 1] = 0
    C.uv[N+1, ] = C.uv[, N+1] = c.uv$x
    names(C.uv) = NULL
    attr(C.uv, "control") <- c(rho = rho)
    C.uv = list(x = U, y = U, z = matrix(C.uv, ncol = length(U)))
    if (!border) {
        C.uv$z = C.uv$z[-1, ]
        C.uv$z = C.uv$z[-N, ]
        C.uv$z = C.uv$z[, -1]
        C.uv$z = C.uv$z[, -N]
        C.uv$x = C.uv$y = V
    }

    # Return Value:
    C.uv
}


# ------------------------------------------------------------------------------


.pellipticalCopulaDiag <-
    function(N, rho = 0.75, param = NULL, type = ellipticalList(), border = TRUE)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Computes elliptical diagonal cross section copula probability

    # Arguments:
    #   see function pellipticalCopula()

    # FUNCTION:

    # Settings:
    U = (0:N)/N
    V = (1:(N-1))/N

    # Compute Density on Regular Grid:
    c.uv = .dellipticalCopulaGrid(N, rho, param, type[1], border = TRUE)
    c.uv$z[is.na(c.uv$z)] = 0

    # Integrate to get Probability:
    C.uu = 0*U
    for (i in 1:(N+1)) {
        C.uu[i] = sum(c.uv$z[1:i, 1:i])
    }
    C.uu = C.uu/N^2
    names(C.uu) = NULL
    attr(C.uu, "control") <- c(rho = rho)

    if (border) {
        C.uu = list(x = U, y = C.uu)
    } else {
        C.uu = list(x = V, y = C.uu[c(-1,-(N+1))])
    }

    # Return Value:
    C.uu
}


# ------------------------------------------------------------------------------


.pellipticalPerspSlider <-
    function(B = 20)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Displays interactively perspective plots of probability

    # Arguments:

    # FUNCTION:

    # Graphic Frame:
    par(mfrow = c(1, 1))

    # Internal Function:
    refresh.code = function(...)
    {
        # Startup Counter:
        .counter <- getRmetricsOptions(".counter") + 1
        setRmetricsOptions(.counter = .counter)
        if (.counter < 7) return ()

        # Sliders:
        Copula = .sliderMenu(no = 1)
        N = .sliderMenu(no = 2)
        rho = .sliderMenu(no = 3)
        nu = .sliderMenu(no = 4)
        s = .sliderMenu(no = 5)
        theta = .sliderMenu(no = 6)
        phi = .sliderMenu(no = 7)
        r = 1

        # Title:
        Names =
            c("- Normal", "- Student t", "- Logistic", "- Exponential Power")
        if (nu == 1) Names[2] = "- Student-t [Cauchy]"
        if (s == 0.5) Names[4] = "- Exponential Power [Laplace]"
        if (s == 1) Names[4] = "- Exponential Power [Kotz|Normal]"
        Title = paste("Elliptical Copula No:", as.character(Copula),
            Names[Copula], "\nrho = ", as.character(rho))
        if (Copula == 2) Title = paste(Title, "nu =", as.character(nu))
        if (Copula == 4) Title = paste(Title, "s =", as.character(s))

        # Plot:
        Type = c("norm", "t", "logistic", "epower")
        param = NULL
        if (Copula == 2) param = nu
        if (Copula == 4) param = c(r, s)
        P = .pellipticalCopulaGrid(N = N, rho = rho, param = param,
            type = Type[Copula], border = TRUE)
        persp(P, theta = theta, phi = phi, col = "steelblue", shade = 0.5,
            ticktype = "detailed", cex = 0.5, xlab = "u", ylab = "v",
            zlab = "C(u, v)", xlim = c(0, 1), ylim = c(0, 1), zlim = c(0, 1) )
        title(main = Title)
        Tau = as.character(round(2*asin(rho)/pi, 2))
        mTitle = paste("Tau", Tau)
        mtext(mTitle, side = 4, col = "grey", cex = 0.7)
        mTitle = paste("1: Normal | 2: Student-t [Cauchy] | 3: Logistic |",
            "4: Exponential Power [Laplace|Kotz]")
        mtext(mTitle, side = 1, line = 3, col = "grey", cex = 0.7)

        # Reset Frame:
        par(mfrow = c(1, 1))
    }

    # Open Slider Menu:
    setRmetricsOptions(.counter = 0)
    plot.names = c("Plot - theta", "... phi")
    .sliderMenu(refresh.code,
        names = c("Copula",  "N", "rho", "2: nu", "4: s", plot.names),
        minima      = c( 1,   10, -0.95,       1,    0.1,  -180,    0),
        maxima      = c( 4,  100,  0.95,       B,      5,   180,  360),
        resolutions = c( 1,   10,  0.05,     0.1,    0.1,     1,    1),
        starts      = c( 1,   20,  0.50,       4,      1,   -40,   30))
}


# ------------------------------------------------------------------------------


.pellipticalContourSlider <- 
    function(B = 20)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Displays interactively perspective plots of probability

    # Arguments:

    # FUNCTION:

    # Graphic Frame:
    par(mfrow = c(1, 1))

    # Internal Function:
    refresh.code =  function(...)
    {
        # Startup Counter:
        .counter <- getRmetricsOptions(".counter") + 1
        setRmetricsOptions(.counter = .counter)
        if (.counter < 7) return ()

        # Sliders:
        Copula = .sliderMenu(no = 1)
        N = .sliderMenu(no = 2)
        rho = .sliderMenu(no = 3)
        nu = .sliderMenu(no = 4)
        s = .sliderMenu(no = 5)
        nlev = .sliderMenu(no = 6)
        ncol = .sliderMenu(no = 7)
        r = 1

        # Title:
        Names =
            c("- Normal", "- Student t", "- Logistic", "- Exponential Power")
        if (nu == 1) Names[2] = "- Student-t [Cauchy]"
        if (s == 0.5) Names[4] = "- Exponential Power [Laplace]"
        if (s == 1) Names[4] = "- Exponential Power [Kotz|Normal]"
        Title = paste("Elliptical Copula No:", as.character(Copula),
            Names[Copula], "\nrho = ", as.character(rho))
        if (Copula == 2) Title = paste(Title, "nu =", as.character(nu))
        if (Copula == 4) Title = paste(Title, "s =", as.character(s))

        # Plot:
        Type = c("norm", "t", "logistic", "epower")
        param = NULL
        if (Copula == 2) param = nu
        if (Copula == 4) param = c(r, s)
        P = .pellipticalCopulaGrid(N = N, rho = rho, param = param,
            type = Type[Copula], border = FALSE)
        image(P, col = heat.colors(ncol), ylab = "v")
        mtext("u", side = 1, line = 2, cex = 0.7)
        contour(P, nlevels = nlev, add = TRUE)
        title(main = Title)
        Tau = as.character(round(2*asin(rho)/pi, 2))
        mTitle = paste("Tau", Tau)
        mtext(mTitle, side = 4, col = "grey", cex = 0.7)
        mTitle = paste("1: Normal | 2: Student-t [Cauchy] | 3: Logistic |",
            "4: Exponential Power [Laplace|Kotz]")
        mtext(mTitle, side = 1, line = 3, col = "grey", cex = 0.7)

        # Reset Frame:
        par(mfrow = c(1, 1))
    }

    # Open Slider Menu:
    setRmetricsOptions(.counter = 0)
    plot.names = c("Plot - levels", "... colors")
    .sliderMenu(refresh.code,
        names = c("Copula", "N", "rho", "2: nu", "4: s", plot.names),
        minima      = c( 1,  10, -0.95,       1,    0.1,     5,   12),
        maxima      = c( 4, 100,  0.95,       B,      5,   100,  256),
        resolutions = c( 1,  10,  0.05,     0.1,    0.1,     5,    4),
        starts      = c( 1,  20,  0.50,       4,      1,    10,   32))
}


################################################################################
# FUNCTION:                  ELLIPTICAL COPULAE DENSITY:
#  dellipticalCopula          Computes elliptical copula density
#  dellipticalSlider          Generates interactive plots of density
#  .dnormCopula               Computes normal copula density
#  .dcauchyCopula             Computes Cauchy copula density
#  .dtCopula                  Computes Student-t copula density
#  .dellipticalCopulaGrid     Fast grid version for elliptical copula density
#  .dellipticalPerspSlider    Interactive perspective plots of density
#  .dellipticalContourSlider  Interactive contour plots of density


dellipticalCopula <-
    function(u = 0.5, v = u, rho = 0.75, param = NULL, type = ellipticalList(),
output = c("vector", "list"), border = TRUE)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Computes extreme value copula density

    # Arguments:
    #   u, v - two numeric values or vectors of the same length at
    #       which the copula will be computed. If 'u' is a list then the
    #       the '$x' and '$y' elements will be used as 'u' and 'v'.
    #       If 'u' is a two column matrix then the first column will
    #       be used as 'u' and the the second as 'v'.
    #   rho - a numeric value setting the coorelation strength, ranging
    #       between minus one and one.
    #   param - additional distributional parameters.
    #   type - the type of the elliptical copula. Either "norm" or
    #       "t" denoting the normal or Student-t copula, respectively.
    #   output - a character string specifying how the output should
    #       be formatted. By default a vector of the same length as
    #       'u' and 'v'. If specified as "list" then 'u' and 'v' are
    #       expected to span a two-dimensional grid as outputted by the
    #       function 'grid2d' and the function returns a list with
    #       elements '$x', 'y', and 'z' which can be directly used
    #       for example by 2D plotting functions.

    # Value:
    #   returns a vector or list of probabilities depending on the
    #   value of the "output" variable.

    # Example:
    #   Diagonal Value: pnormCopula((0:10)/10)
    #   persp(pnormCopula(u = grid2d(), output = "list"))

    # FUNCTION:

    # Use Grid Version?
    if (is.numeric(u)) {
        if (length(u) == 1 & u[1] > 1) {
            ans = .dellipticalCopulaGrid(N = u, rho = rho, param = param,
                type = type, border = border)
            return(ans)
        }
    }

    # Match Arguments:
    type = match.arg(type)
    output = match.arg(output)

    # Settings:
    if (is.list(u)) {
        v = u[[2]]
        u = u[[1]]
    }
    if (is.matrix(u)) {
        v = u[, 2]
        u = u[, 1]
    }
    if (length(u) == 1 & u[1] > 1) {
        return(.pellipticalCopulaGrid(N = u, rho, param, type, border = border))
    }

    # Parameters:
    if (type == "t") if (is.null(param)) param = c(nu = 4)
    if (type == "kotz") if (is.null(param)) param = c(r = 1)
    if (type == "epower") if (is.null(param)) param = c(r = 1, s = 1)

    # Density:
    x = .qelliptical(u, param = param, type = type)
    y = .qelliptical(v, param = param, type = type)
    c.uv = delliptical2d(x, y, rho = rho, param = param, type = type) / (
        .delliptical(x, param = param, type = type) *
        .delliptical(y, param = param, type = type) )
    if (rho == 0 & type == "norm") c.uv[!is.na(c.uv)] = 1
    names(c.uv) = NULL
    attr(c.uv, "control") <- c(rho = rho)
    if (output == "list") {
        N = sqrt(length(u))
        x = u[1:N]
        y = matrix(v, ncol = N)[1, ]
        c.uv = list(x = x, y = y, z = matrix(c.uv, ncol = N))
    }

    # Return Value:
    c.uv
}


# ------------------------------------------------------------------------------


dellipticalSlider <- 
    function(type = c("persp", "contour"), B = 20)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Displays interactively plots of density

    # Description:
    #   Displays interactively plots of density

    # Arguments:
    #   type - a character string specifying the plot type.
    #       Either a perspective plot which is the default or
    #       a contour plot with an underlying image plot will
    #       be created.
    #   B - the maximum slider menu value when the boundary
    #       value is infinite. By default this is set to 10.

    # FUNCTION:

    # Settings:
    type = match.arg(type)

    # Plot:
    if (type == "persp")
        .dellipticalPerspSlider(B = B)
    if (type == "contour")
        .dellipticalContourSlider(B = B)

    # Return Value:
    invisible()
}


# ------------------------------------------------------------------------------


.dnormCopula <-
    function(u = 0.5, v = u, rho = 0.75, output = c("vector", "list") )
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Computes normal copula density

    # Arguments:
    #   see function 'dellipticalCopula'

    # FUNCTION:

    # Type:
    output = match.arg(output)

    # Settings:
    type = "norm"
    if (is.list(u)) {
        v = u[[2]]
        u = u[[1]]
    }
    if (is.matrix(u)) {
        v = u[, 2]
        u = u[, 1]
    }

    # Copula Density:
    x = qnorm(u)
    y = qnorm(v)
    c.uv = dnorm2d(x, y, rho)/(dnorm(x) * dnorm(y))
    names(c.uv) = NULL

    # Result:
    attr(c.uv, "control") <- c(rho = rho)

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


.dtCopula <- 
    function(u = 0.5, v = u, rho = 0.75, nu = 4, output = c("vector", "list") )
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Computes Student-t copula density

    # Arguments:
    #   see function 'dellipticalCopula'

    # FUNCTION:

    # Match Arguments:
    output = match.arg(output)

    # Settings:
    type = "t"
    if (is.list(u)) {
        v = u[[2]]
        u = u[[1]]
    }
    if (is.matrix(u)) {
        v = u[, 2]
        u = u[, 1]
    }

    # Copula Probability:
    x = qt(u, df = nu)
    y = qt(v, df = nu)
    c.uv = dt2d(x, y, rho, nu)/(dt(x, nu) * dt(y, nu))
    names(c.uv) = NULL

    # Result:
    attr(c.uv, "control") <- c(rho = rho, nu = nu)

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


.dcauchyCopula <- 
    function(u = 0.5, v = u, rho = 0.75, nu = 4, output = c("vector", "list") )
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Computes Student-t copula density

    # Arguments:
    #   see function 'dellipticalCopula'

    # FUNCTION:

    # Cauchy Density:
    c.uv = .dtCopula(u = u, v = v, rho = rho, nu = 1, output = output)
    attr(c.uv, "control") <- c(rho = rho)

    # Return Value:
    c.uv
}


# ------------------------------------------------------------------------------


.dellipticalCopulaGrid <- 
    function(N, rho = 0.75, param = NULL, type = ellipticalList(), border = TRUE)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Computes extreme value copula density

    # Arguments:
    #   N - the number of grid points is (N+1)*(N+1)
    #   rho - a numeric value setting the coorelation strength, ranging
    #       between minus one and one.
    #   param - additional distributional parameters.
    #   type - the type of the elliptical copula. Either "norm" or
    #       "t" denoting the normal or Student-t copula, respectively.

    # Value:
    #   returns a vector or list of probabilities depending on the
    #   value of the "output" variable.

    # Note:
    #   Made for the Sliders.

    # FUNCTION:

    # Settings:
    type = type[1]
    U = (0:N)/N
    V = (1:(N-1))/N

    # Reduce to Grid - speeds up the computation:
    M = N%/%2 + 1
    X = .qelliptical(U[1:M], param = param, type = type)
    if (N%%2 == 0) {
        X = c(X, rev(-X)[-1])
    } else {
        X = c(X, rev(-X))
    }
    NX = length(X)
    x = rep(X, times = NX)
    y = rep(X, each = NX)
    D = .delliptical(X, param = param, type = type)
    DX = rep(D, times = NX)
    DY = rep(D, each = NX)

    # Density:
    c.uv = delliptical2d(x, y, rho = rho, param = param, type = type) / (DX*DY)
    if (rho == 0 & type == "norm") c.uv[!is.na(c.uv)] = 1
    c.uv[is.na(c.uv)] = 0
    names(c.uv) = NULL
    attr(c.uv, "control") <- c(rho = rho)
    c.uv = list(x = U, y = U, z = matrix(c.uv, ncol = N+1))
    if (!border) {
        c.uv$z = c.uv$z[-1, ]
        c.uv$z = c.uv$z[-N, ]
        c.uv$z = c.uv$z[, -1]
        c.uv$z = c.uv$z[, -N]
        c.uv = list(x = V, y = V, z = matrix(c.uv$z, ncol = N-1))
    }

    # Return Value:
    c.uv
}


# ------------------------------------------------------------------------------


.dellipticalPerspSlider <- 
    function(B = 20)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Displays interactively perspective plots of density

    # FUNCTION:

    # Graphic Frame:
    par(mfrow = c(1, 1))

    # Internal Function:
    refresh.code =  function(...)
    {
        # Startup Counter:
        .counter <- getRmetricsOptions(".counter") + 1
        setRmetricsOptions(.counter = .counter)
        if (.counter < 7) return ()

        # Sliders:
        Copula = .sliderMenu(no = 1)
        N = .sliderMenu(no = 2)
        rho = .sliderMenu(no = 3)
        nu = .sliderMenu(no = 4)
        s = .sliderMenu(no = 5)
        theta = .sliderMenu(no = 6)
        phi = .sliderMenu(no = 7)
        r = 1

        # Title:
        Names =
            c("- Normal", "- Student t", "- Logistic", "- Exponential Power")
        if (nu == 1) Names[2] = "- Student-t [Cauchy]"
        if (s == 0.5) Names[4] = "- Exponential Power [Laplace]"
        if (s == 1) Names[4] = "- Exponential Power [Kotz|Normal]"
        Title = paste("Elliptical Copula Density No:", as.character(Copula),
            Names[Copula], "\nrho = ", as.character(rho))
        if (Copula == 2) Title = paste(Title, "nu =", as.character(nu))
        if (Copula == 4) Title = paste(Title, "s =", as.character(s))

        # Plot:
        uv = grid2d(x = (1:(N-1))/N)
        Type = c("norm", "t", "logistic", "epower")
        param = NULL
        if (Copula == 2) param = nu
        if (Copula == 4) param = c(r, s)
        D = .dellipticalCopulaGrid(N, rho = rho, param = param,
            type = Type[Copula], border = FALSE)
        Integrated = as.character(round(mean(D$z),2))
        Var = var(as.vector(D$z), na.rm = TRUE)
        if (Var < 1.0e-6) {
            # A flat perspective plot fails, if zlim is not specified!
            Mean = round(1.5*mean(as.vector(D$z), na.rm = TRUE), 2)
            persp(D, theta = theta, phi = phi, col = "steelblue", shade = 0.5,
                ticktype = "detailed", cex = 0.5, xlab = "u", ylab = "v",
                zlim = c(0, Mean), zlab = "C(u,v)" )
        } else {
            persp(D, theta = theta, phi = phi, col = "steelblue", shade = 0.5,
                ticktype = "detailed", cex = 0.5, xlab = "u", ylab = "v",
                zlab = "C(u,v)" )
        }
        title(main = Title)
        Tau = as.character(round(2*asin(rho)/pi, 2))
        mTitle = paste("Mean: ", Integrated, " |  Tau", Tau)
        mtext(mTitle, side = 4, col = "grey", cex = 0.7)
        mTitle = paste("1: Normal | 2: Student-t [Cauchy] | 3: Logistic |",
            "4: Exponential Power [Laplace|Kotz]")
        mtext(mTitle, side = 1, col = "grey", cex = 0.7)

        # Reset Frame:
        par(mfrow = c(1, 1))
    }

    # Open Slider Menu:
    setRmetricsOptions(.counter = 0)
    plot.names = c("Plot - theta", "... phi")
    .sliderMenu(refresh.code,
        names = c("Copula", "N", "rho", "3: nu", "4: s", plot.names),
        minima      = c( 1,  10, -0.95,       1,   0.1,  -180,    0),
        maxima      = c( 4, 100,  0.95,       B,     5,   180,  360),
        resolutions = c( 1,  10,  0.05,     0.1,   0.1,     1,    1),
        starts      = c( 1,  20,  0.50,       4,     1,   -40,   30))
}


# ------------------------------------------------------------------------------


.dellipticalContourSlider <-
    function(B = 20)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Displays interactively perspective plots of density

    #FUNCTION:

    # Graphic Frame:
    par(mfrow = c(1, 1))

    # Internal Function:
    refresh.code =  function(...)
    {
        # Startup Counter:
        .counter <- getRmetricsOptions(".counter") + 1
        setRmetricsOptions(.counter = .counter)
        if (.counter < 7) return ()

        # Sliders:
        Copula = .sliderMenu(no = 1)
        N = .sliderMenu(no = 2)
        rho = .sliderMenu(no = 3)
        nu = .sliderMenu(no = 4)
        s = .sliderMenu(no = 5)
        nlev = .sliderMenu(no = 6)
        ncol = .sliderMenu(no = 7)
        if (rho == 0 & Copula == 1) return(invisible())
        r = 1

        # Title:
        Names =
            c("- Normal", "- Student t", "- Logistic", "- Exponential Power")
        if (nu == 1) Names[2] = "- Student-t [Cauchy]"
        if (s == 0.5) Names[4] = "- Exponential Power [Laplace]"
        if (s == 1) Names[4] = "- Exponential Power [Kotz|Normal]"
        Title = paste("Elliptical Copula Density No:", as.character(Copula),
            Names[Copula], "\nrho = ", as.character(rho))
        if (Copula == 2) Title = paste(Title, "nu =", as.character(nu))
        if (Copula == 4) Title = paste(Title, "s =", as.character(s))

        # Plot:
        uv = grid2d(x = (0:N)/N)
        Type = c("norm", "t", "logistic", "laplace", "kotz", "epower")
        param = NULL
        if (Copula == 2) param = nu
        if (Copula == 5) param = r
        if (Copula == 6) param = c(r, s)
        D = .dellipticalCopulaGrid(N, rho = rho, param = param,
            type = Type[Copula], border = FALSE)
        Integrated = as.character(round(mean(D$z),2))
        image(D, col = heat.colors(ncol), ylab = "v",
            xlim = c(0,1), ylim = c(0,1) )
        mtext("u", side = 1, line = 2, cex = 0.7)
        contour(D, nlevels = nlev, add = TRUE)
        title(main = Title)
        Tau = as.character(round(2*asin(rho)/pi, 2))
        mTitle = paste("Mean: ", Integrated, " |  Tau", Tau)
        mtext(mTitle, side = 4, col = "grey", cex = 0.7)
        mTitle = paste("1: Normal | 2: Student-t [Cauchy] | 3: Logistic |",
            "4: Exponential Power [Laplace|Kotz]")
        mtext(mTitle, side = 1, line =  3, col = "grey", cex = 0.7)

        # Reset Frame:
        par(mfrow = c(1, 1))
    }

    # Open Slider Menu:
    setRmetricsOptions(.counter = 0)
    plot.names = c("Plot - levels", "... colors")
    .sliderMenu(refresh.code,
        names = c("Copula", "N", "rho", "2: nu", "4: s", plot.names),
        minima      = c( 1,  10, -0.95,       1,    0.1,     5,   12),
        maxima      = c( 4, 100,  0.95,       B,      5,   100,  256),
        resolutions = c( 1,  10,  0.05,     0.1,    0.1,     5,    4),
        starts      = c( 1,  20,  0.50,       4,      1,    10,   32))
}


################################################################################

