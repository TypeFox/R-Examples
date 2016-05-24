
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

# Copyrights (C)
# for this R-port:
#   1999 - 2009, Diethelm Wuertz, GPL
#   Diethelm Wuertz <wuertz@itp.phys.ethz.ch>
#   info@rmetrics.org
#   www.rmetrics.org
# for the code accessed (or partly included) from other R-ports:
#   see R's copyright and license files
# for the code accessed (or partly included) from contributed R-ports
# and other sources
#   see Rmetrics's copyright file


################################################################################
# FUNCTION:             GEV DISTRIBUTION FAMILY: [CALLING EVD]
#  dgev                  Density for the GEV Distribution
#   pgev                  Probability for the GEV Distribution
#   qgev                  Quantiles for the GEV Distribution
#   rgev                  Random variates for the GEV Distribution
#  gevMoments            Computes true statistics for GEV distribution
#  gevSlider             Displays distribution and rvs for GEV distribution
# FUNCTION:             GEV DISTRIBUTION FAMILY: [USE FROM EVD]
#  .devd                 Density for the GEV Distribution
#   .pevd                 Probability for the GEV Distribution
#   .qevd                 Quantiles for the GEV Distribution
#   .revd                 Random variates for the GEV Distribution
################################################################################


dgev =
function(x, xi = 1, mu = 0, beta = 1, log = FALSE)
{   
    # A function implemented from package evd

    # Description:
    #   GEV Density Function
    #   Note: 1 + xi*(x-mu)/beta > 0
    #   xi > 0 Frechet
    #   xi = 0 Gumbel
    #   xi < 0 weibl

    # FUNCTION:

    # Density:
    d = .devd(x, location = mu, scale = beta, shape = xi, log = FALSE)

    # Add Attribute:
    attr(d, "control") = data.frame(xi = xi, mu = mu, beta = beta,
        log = log, row.names = "")

    # Return Value:
    d
}


# ------------------------------------------------------------------------------


pgev =
function(q, xi = 1, mu = 0, beta = 1, lower.tail = TRUE)
{   
    # A function implemented from package evd

    # Description:
    #   GEV Probability Function
    #   Note: 1 + xi*(x-mu)/beta > 0
    #   xi > 0 Frechet
    #   xi = 0 Gumbel
    #   xi < 0 Weibull

    # FUNCTION:

    # Probability:
    p = .pevd(q, location = mu, scale = beta, shape = xi,
        lower.tail = lower.tail)

    # Add Attribute:
    attr(p, "control") = data.frame(xi = xi, mu = mu, beta = beta,
        lower.tail = lower.tail, row.names = "")

    # Return Value:
    p
}


# ------------------------------------------------------------------------------


qgev =
function(p, xi = 1, mu = 0, beta = 1, lower.tail = TRUE)
{   
    # A function implemented from package evd

    # Description:
    #   GEV Quantile Function
    #   Note: 1 + xi*(x-mu)/beta > 0
    #   xi > 0 Frechet
    #   xi = 0 Gumbel
    #   xi < 0 Weibull

    # FUNCTION:

    # Quantiles:
    q = .qevd(p, location = mu, scale = beta, shape = xi,
        lower.tail = lower.tail)

    # Add Attribute:
    attr(q, "control") = data.frame(xi = xi, mu = mu, beta = beta,
        lower.tail = lower.tail, row.names = "")

    # Return Value:
    q
}


# ------------------------------------------------------------------------------


rgev =
function(n, xi = 1, mu = 0, beta = 1)
{   
    # A function implemented from package evd

    # Description:
    #   GEV Random Variables
    #   Note: 1 + xi*(x-mu)/beta > 0
    #   xi > 0 Frechet
    #   xi = 0 Gumbel
    #   xi < 0 Weibull

    # FUNCTION:

    # Random Variates:
    r = .revd(n = n, location = mu, scale = beta, shape = xi)

    # Add Attribute:
    attr(q, "control") = data.frame(xi = xi, mu = mu, beta = beta)

    # Return Value:
    r
}


# ----------------------------------------------------------------------------


gevMoments =
function(xi = 0, mu = 0, beta = 1)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Compute true statistics for Generalized Extreme Value distribution

    # Value:
    #   Returns true mean for xi < 1 and variance for xi < 1/2
    #   of GEV distribution, otherwise NaN is returned

    # FUNCTION:

    # MEAN: Returns for x >= 1 NaN:
    g = c(1, 0, NaN)
    xinv = 1/ ( xi + sign(abs(xi)) - 1 )

    # For xi = the result is eulers constant
    euler = 0.57721566490153286060651209008240243104
    xi0 = c(0, beta*euler, 0)

    # Supress warning for NaN's from Gamma Function:
    warn.old <- getOption("warn")
    options(warn = -1)
    gevMean = mu + beta * xinv * (gamma(1-xi)-1) * g[sign(xi-1)+2] +
        xi0[(sign(xi)+2)]
    options(warn = warn.old)

    # VAR: Returns for x >= 1 NaN:
    g = c(1, 0, NaN)
    xinv = 1/ ( xi + sign(abs(xi)) - 1 )
    xi0 = c(0, (beta*pi)^2 / 6, 0)

    # Supress warning for NaN's from Gamma Function:
    warn.old <- getOption("warn")
    options(warn = -1)
    gevVar = (beta*xinv)^2 * (gamma(1-2*xi) - gamma(1-xi)^2 ) *
        g[sign(2*xi-1)+2] + xi0[(sign(xi)+2)]
    options(warn = warn.old)

    # Result:
    param = c(xi = xi, mu = mu, beta = beta)
    ans = list(param = param, mean = gevMean, var = gevVar)

    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


gevSlider =
function(method = c("dist", "rvs"))
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Displays distribution and rvs for GEV distribution

    # FUNCTION:

    # Settings:
    method = match.arg(method)

    # Internal Function:
    refresh.code = function(...)
    {
        # Sliders:
        N      = .sliderMenu(no = 1)
        xi     = .sliderMenu(no = 2)
        mu     = .sliderMenu(no = 3)
        beta   = .sliderMenu(no = 4)

        # Compute Data:
        if (xi > 0) {
            pmin = 0
            pmax = 0.999
        } else if (xi < 0) {
            pmin = 0.001
            pmax = 1
        } else if (xi == 0) {
            pmin = 0.001
            pmax = 0.999
        }
        xmin = round(qgev(pmin, xi, mu, beta), digits = 2)
        xmax = round(qgev(pmax, xi, mu, beta), digits = 2)
        s = seq(xmin, xmax, length = N)
        y1 = dgev(s, xi, mu, beta)
        y2 = pgev(s, xi, mu, beta)
        Moments = gevMoments(xi, mu, beta)
        Mean = round(Moments$mean, 2)
        Var = round(Moments$var, 2)
        mText = paste("Mean =", Mean, " | Variance = ", Var)
        main1 = paste("GEV Density\n",
            "xi = ", as.character(xi), " | ",
            "mu = ", as.character(mu), " | ",
            "beta = ", as.character(beta) )
        main2 = paste("GEV Probability\n",
            "xmin [0.00] = ", as.character(xmin), " | ",
            "xmax [0.99] = ", as.character(xmax) )
        Median = qgev(0.5, xi, mu, beta)

        # Frame:
        par(mfrow = c(2, 1), cex = 0.7)

        # Density:
        if (method == "rvs") {
            x = rgev(N, xi, mu, beta)
            hist(x, probability = TRUE, col = "steelblue", border = "white",
                xlim = c(xmin, xmax), ylim = c(0, 1.1*max(y1)), main = main1,
                breaks = "FD" )
            lines(s, y1, col = "orange")
            mtext(mText, side = 4, col = "grey", cex = 0.7)
        } else {
            plot(s, y1, type = "l", xlim = c(xmin, xmax), col = "steelblue")
            abline(h = 0, lty = 3)
            abline(v = Median, lty = 3, col = "red")
            abline(v = Mean, lty = 3, col = "darkgreen")
            title(main = main1)
            mtext(mText, side = 4, col = "grey", cex = 0.7)
        }

        # Probability:
        plot(s, y2, type = "l", xlim = c(xmin, xmax), ylim = c(0, 1),
            col = "steelblue" )
        abline(h = 0, lty = 3)
        abline(h = 0.5, lty = 3, col = "red")
        abline(v = Median, lty = 3, col = "red")
        abline(v = Mean, lty = 3, col = "darkgreen")
        title(main = main2)
        mtext(mText, side = 4, col = "grey", cex = 0.7)

        # Reset Frame:
        par(mfrow = c(1, 1), cex = 0.7)
    }

    # Open Slider Menu:
    .sliderMenu(refresh.code,
       names =       c(   "N", "xi",  "mu",  "beta"),
       minima =      c(   50, -1.50, -5.00,   0.10 ),
       maxima =      c( 1000,  1.50, +5.00,   5.00 ),
       resolutions = c(   50,  0.01,  0.10,   0.10 ),
       starts =      c(  500,  0.00,  0.00,   1.00 )
    )
}



################################################################################


.devd =
function(x, location = 0, scale = 1, shape = 0, log = FALSE)
{   # A modified copy from contributed R package evd

    # FUNCTION:

    # Check:
    stopifnot(min(scale) > 0)
    stopifnot(length(shape) == 1)

    # Density:
    x = (x - location) / scale
    if (shape == 0) {
        d = log(1/scale) - x - exp(-x)
    } else {
        nn = length(x)
        xx = 1 + shape * x
        xxpos = xx[xx > 0 | is.na(xx)]
        scale = rep(scale, length.out = nn)[xx > 0 | is.na(xx)]
        d = numeric(nn)
        d[xx > 0 | is.na(xx)] = log(1/scale) - xxpos^(-1/shape) -
            (1/shape + 1) * log(xxpos)
        d[xx <= 0 & !is.na(xx)] = -Inf
    }

    # Log:
    if (!log) {
        d = exp(d)
    }

    # Add Attribute:
    attr(d, "control") = data.frame(location = location[1], scale = scale[1],
        shape = shape[1], log = log, row.names = "")


    # Return Value:
    d
}


# ------------------------------------------------------------------------------


.pevd =
function(q, location = 0, scale = 1, shape = 0, lower.tail = TRUE)
{   # A modified copy from contributed R package evd

    # FUNCTION:

    # Check:
    stopifnot(min(scale) > 0)
    stopifnot(length(shape) == 1)

    # Probabilities:
    q = (q - location)/scale
    if (shape == 0) {
        p = exp(-exp(-q))
    } else {
        p = exp(-pmax(1 + shape * q, 0)^(-1/shape))
    }

    # Lower Tail:
    if (!lower.tail) {
        p = 1 - p
    }

    # Add Attribute:
    attr(p, "control") = data.frame(location = location[1], scale = scale[1],
        shape = shape[1], lower.tail = lower.tail, row.names = "")

    # Return Value:
    p
}


# ------------------------------------------------------------------------------


.qevd =
function(p, location = 0, scale = 1, shape = 0, lower.tail = TRUE)
{   # A modified copy from contributed R package evd

    # FUNCTION:

    # Check:
    stopifnot(min(scale) > 0)
    stopifnot(length(shape) == 1)
    stopifnot(min(p, na.rm = TRUE) >= 0)
    stopifnot(max(p, na.rm = TRUE) <= 1)

    # Quantiles:
    if (!lower.tail) p = 1 - p
    if (shape == 0) {
        q = location - scale * log(-log(p))
    } else {
        q = location + scale * ((-log(p))^(-shape) - 1)/shape
    }

    # Add Attribute:
    attr(q, "control") = data.frame(location = location[1], scale = scale[1],
        shape = shape[1], lower.tail, row.names = "")

    # Return Value:
    q
}


# ------------------------------------------------------------------------------


.revd =
function(n, location = 0, scale = 1, shape = 0)
{   # A modified copy from contributed R package evd

    # FUNCTION:

    # Check:
    stopifnot(min(scale) > 0)
    stopifnot(length(shape) == 1)

    # Random Variates:
    if (shape == 0) {
        r = location - scale * log(rexp(n))
    } else {
        r = location + scale * (rexp(n)^(-shape) - 1)/shape
    }

    # Add Attribute:
    attr(r, "control") = data.frame(location = location[1], scale = scale[1],
        shape = shape[1], row.names = "")

    # Return Value:
    r
}


################################################################################

