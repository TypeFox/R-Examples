
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA


################################################################################
# FUNCTION:             BIVARIATE NORMAL DISTRIBUTION:
#  dnorm2d               Computes bivariate Normal density function
#  pnorm2d               Computes bivariate Normal probability function
#  rnorm2d               Generates bivariate normal random deviates
################################################################################


pnorm2d <-
    function(x, y = x, rho = 0)
{   
    # pnorm2d: A copy from R package "sn"

    # Description:
    #   Computes bivariate Normal probability function

    # Arguments:
    #   x, y - two numeric values or vectors of the same length at
    #       which the probability will be computed.

    # Value:
    #   returns a numeric vector of probabilities of the same length
    #   as the input vectors

    # FUNCTION:

    # Probaility:
    X <- cbind(x, y)
    ans <- apply(X, 1, .pnorm2d, rho = rho)
    attr(ans, "control") = c(rho = rho)


    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


.pnorm2d <- 
    function(X, rho = 0)
{   
    # pnorm2d: A copy from R package "sn"

    # Description:
    #   Bivariate Normal probability function

    # Arguments:
    #   x, y - two numeric values at which the probability will
    #   be computed.

    # Value:
    #   returns a numeric vector of probabilities of the same length
    #   as the input vectors

    # FUNCTION:

    # Probability:
    x <- X[1]
    y <- X[2]
    if (x == 0 & y == 0) {
        return(0.25 + asin(rho)/(2 * pi))
    }
    p = 0.5 * (pnorm(x) + pnorm(y))
    if (x == 0) {
        p = p - 0.25 * sign(y)
    } else {
        if (is.finite(x)) {
            Y = (y - rho * x)/(x * sqrt(1 - rho^2))
        } else {
            Y = -rho/sqrt(1-rho^2)
        }
        p = p - .TOwen(x, Y)
    }
    if (y == 0) {
        p = p - 0.25 * sign(x)
    } else {
        if (is.finite(y)) {
            X = (x - rho * y)/(y * sqrt(1 - rho^2))
        } else {
            X = -rho/sqrt(1-rho^2)
        }
        p = p - .TOwen(y, X)
    }
    if (is.finite(x) & is.finite(y)) {
        if ((x * y < 0) | ((x * y == 0) & (x + y) < 0)) {
            p = p - 0.5
        }
    }

    # Return Value:
    return(p)
}


# ------------------------------------------------------------------------------


.TInt =
    function(h, a, jmax, cut.point)
{   
    # T.int: A copy from R package "sn"

    # Note:
    #   Required by .pnorm2d and .TOwen

    # FUNCTION:

    .fui = function(h, i) (h^(2 * i))/((2^i) * gamma(i + 1))
    seriesL = seriesH = NULL
    i = 0:jmax
    low = (h <= cut.point)
    hL = h[low]
    hH = h[!low]
    L = length(hL)
    if (L > 0) {
        b = outer(hL, i, .fui)
        cumb = apply(b, 1, cumsum)
        b1 = exp(-0.5 * hL^2) * t(cumb)
        matr = matrix(1, jmax + 1, L) - t(b1)
        jk = rep(c(1, -1), jmax)[1:(jmax + 1)]/(2 * i + 1)
        matr = t(matr * jk) %*% a^(2 * i + 1)
        seriesL = (atan(a) - as.vector(matr))/(2 * pi)
    }
    if (length(hH) > 0) {
        seriesH = atan(a) * exp(-0.5 * (hH^2) * a/atan(a)) *
            (1 + 0.00868 * (hH^4) * a^4)/(2 * pi)
    }
    series = c(seriesL, seriesH)
    id = c((1:length(h))[low], (1:length(h))[!low])
    series[id] = series

    # Return Value:
    series
}


# ------------------------------------------------------------------------------


.TOwen <- 
    function (h, a, jmax = 50, cut.point = 6)
{   
    # T.Owen: A copy from R package "sn"

    # Note:
    #   Required by .pnorm2d

    # FUNCTION:

    if (!is.vector(a) | length(a) > 1)
        stop("a must be a vector of length 1")
    if (!is.vector(h))
        stop("h must be a vector")
    aa = abs(a)
    ah = abs(h)
    if (aa == Inf)
        return(0.5 * pnorm(-ah))
    if (aa == 0)
        return(rep(0, length(h)))
    na = is.na(h)
    inf = (ah == Inf)
    ah = replace(ah, (na | inf), 0)
    if (aa <= 1) {
        owen = .TInt(ah, aa, jmax, cut.point)
    } else {
        owen = 0.5 * pnorm(ah) + pnorm(aa * ah) * (0.5 - pnorm(ah)) -
            .TInt(aa * ah, (1/aa), jmax, cut.point)
    }
    owen = replace(owen, na, NA)
    owen = replace(owen, inf, 0)
    ans = return(owen * sign(a))

    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


dnorm2d <- 
    function(x, y = x, rho = 0)
{   
    # A function implemented by Diethelm Wuertz

    # Arguments:
    #   x,y - two numeric vectors
    #   rho - the linear correlation, a numeric value between
    #       minus one and one.

    # FUNCTION:

    # Argument:
    xoy = (x^2 - 2*rho*x*y + y^2)/ (2*(1 - rho^2))

    # Density:
    density = exp(-xoy) / ( 2*pi*sqrt(1-rho^2))
    attr(density, "control") = c(rho = rho)

    # Return Value:
    density
}


# ------------------------------------------------------------------------------


.dnorm2d =
function(x, y = x, rho = 0)
{   # A function implemented by Diethelm Wuertz

    # Arguments:
    #   x,y - two numeric vectors
    #   rho - the linear correlation, a numeric value between
    #       minus one and one.

    # Note:
    #   Partly copied from contributed R package 'mvtnorm'
    #   Author Friedrich Leisch

    # FUNCTION

    # Settings:
    mean = c(0,0)
    sigma = diag(2)
    sigma[1,2] = sigma[2,1] = rho
    log = FALSE
    x = cbind(x, y)

    # From mvtnorm - Check:
    if (is.vector(x)) {
        x = matrix(x, ncol = length(x))
    }
    if (missing(mean)) {
        mean = rep(0, length = ncol(x))
    }
    if (missing(sigma)) {
        sigma = diag(ncol(x))
    }
    if (ncol(x) != ncol(sigma)) {
        stop("x and sigma have non-conforming size")
    }
    if (nrow(sigma) != ncol(sigma)) {
        stop("sigma meanst be a square matrix")
    }
    if (length(mean) != nrow(sigma)) {
        stop("mean and sigma have non-conforming size")
    }

    # From mvtnorm - Density:
    distval = mahalanobis(x, center = mean, cov = sigma)
    logdet = sum(log(eigen(sigma, symmetric = TRUE, only.values = TRUE)$values))
    logretval = -(ncol(x)*log(2*pi) + logdet + distval)/2
    if(log) return(logretval)
    ans = exp(logretval)
    attr(ans, "control") = c(rho = rho)

    # Return value:
    ans
}


# ------------------------------------------------------------------------------


rnorm2d <- 
    function(n, rho = 0)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Generates bivariate normal random deviates

    # Arguments:
    #   n - number of random deviates to be generated
    #   rho - the linear correlation, a numeric value between
    #       minus one and one.

    # Note:
    #   Partly copied from contributed R package 'mvtnorm'
    #   Author Friedrich Leisch

    # FUNCTION

    # Settings:
    mean = c(0,0)
    sigma = diag(2)
    sigma[1,2] = sigma[2,1] = rho

    # From mvtnorm - Random Numbers:
    ev = eigen(sigma, symmetric = TRUE)$values
    if (!all(ev >= -sqrt(.Machine$double.eps) * abs(ev[1])))
        warning("sigma is numerically not positive definite")
    sigsvd = svd(sigma)
    ans = t(sigsvd$v %*% (t(sigsvd$u) * sqrt(sigsvd$d)))
    ans = matrix(rnorm(n * ncol(sigma)), nrow = n) %*% ans
    ans = sweep(ans, 2, mean, "+")
    attr(ans, "control") = c(rho = rho)

    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


.rnorm2d <- 
    function(n, rho = 0)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Alternative direct algorithm from Lindskog Master Thesis

    # Arguments:
    #   n - number of random deviates to be generated
    #   rho - the linear correlation, a numeric value between
    #       minus one and one.

    # FUNCTION:

    # Random Deviates
    x = matrix(c(1, rho, rho,1), 2)
    V = NULL
    U = chol(x)
    siz = dim(x)[1]
    for(i in 1:n) {
        Z = rnorm(siz)
        res = t(U)%*%Z
        V = cbind(V,res)
    }
    rmn = t(V)

    # Return Value:
    rmn
}


###############################################################################


