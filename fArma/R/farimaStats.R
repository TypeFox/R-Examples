
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# A copy of the GNU General Public License is available via WWW at
# http://www.gnu.org/copyleft/gpl.html.  You can also obtain it by
# writing to the Free Software Foundation, Inc., 59 Temple Place,
# Suite 330, Boston, MA  02111-1307  USA.

# Copyrights (C)
# for this R-port:
#   1999 - 2007, Diethelm Wuertz, GPL
#   Diethelm Wuertz <wuertz@itp.phys.ethz.ch>
#   info@rmetrics.org
#   www.rmetrics.org
# for the code accessed (or partly included) from other R-ports:
#   see R's copyright and license files
# for the code accessed (or partly included) from contributed R-ports
# and other sources
#   see Rmetrics's copyright file


################################################################################
# FUNCTIONS:            DESCRIPTION:
#  farimaTrueacf         Returns FARMA true autocorrelation function
#  farimaTruefft         Returns FARMA true fast Fourier transform
#  .farimaStatsSlider     Displays farima true statistics
################################################################################


################################################################################
# DESCRIPTION:
#   The functions are from the appendix of J. Beran "Statistics for
#   long-memory processes", Chapman and Hall 1984
# LICENSE:
#   Reimplemented functions from Beran's SPlus Scripts.
#   Permission is hereby given to StatLib to redistribute this software.
#   The software can be freely used for non-commercial purposes, and can
#   be freely distributed for non-commercial purposes only.
# AUTHORS:
#   Jan Beran <jberan@iris.rz.uni-konstanz.de>
#   Modified: Martin Maechler <maechler@stat.math.ethz.ch>
#   Modified: Diethelm Wuertz <wuertz@itp.phys.ethz.ch> for this R-Port


# ------------------------------------------------------------------------------


farimaTrueacf =
function(n = 100, H = 0.7)
{   # A function implemented by Diethelm Wuertz

    # FUNCTION:

    # ACF:
    ans = .ckFARIMA0(n = n, H = H)

    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


farimaTruefft =
function(n = 100, H = 0.7)
{   # A function implemented by Diethelm Wuertz

    # FUNCTION:

    # FFT:
    ans = .gkFARIMA0(n = n, H = H)

    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


.ckFARIMA0 =
function(n, H)
{
    # Description:
    #   Computes the covariances of a fractional ARIMA(0,d,0) process

    # Arguments:
    #   n = length of time series
    #   H = self-similarity parameter

    # Value:
    #   Covariances up to lag n-1

    # Author:
    #   Jan Beran; modified: Martin Maechler, Date: Sep 95.

    # FUNCTION:

    # Covariances:
    # result = (0:(n-1))
    # k = 1:(n-1)
    # d = H - 0.5
    # result[1] = gamma(1-2*d)/gamma(1-d)**2
    # result[k+1] = result[1]*(k**(2*H-2))*gamma(1-d)/gamma(d)
    # result[k+1] = result[1]*gamma(k+d)* gamma(1-d)/(gamma(k-d+1)*gamma(d))
    # ans = drop(result)

    # Covariances:
    res = numeric(n)
    d = H - 0.5
    g1d = gamma(1 - d)
    gd = pi/(sin(pi * d) * g1d)
    res[1] = gamma(1 - 2 * d)/g1d^2
    k = 1:min(50, n - 1)
    res[k + 1] = res[1] * gamma(k + d) * g1d/(gamma(k - d + 1) * gd)
    if (n > 51) {
        k <- 51:(n - 1)
        res[k + 1] <- res[1] * g1d/gd * k^(2 * H - 2)
    }

    # Return Value:
    res
}


# ------------------------------------------------------------------------------


farimaStatsSlider =
function()
{   # A function implemented by Diethelm Wuertz

    # Description
    #   Displays farima true Statistics: ACF and FFT

    # Example:
    #   farimaStatsSlider()

    # FUNCTION:

    # Internal Function:
    refresh.code = function(...)
    {
        # Sliders:
        n = fBasics:::.sliderMenu(no = 1)
        H = fBasics:::.sliderMenu(no = 2)

        # Frame:
        par(mfrow = c(2, 1), cex = 0.7)

        # FGN ACF:
        ans = farimaTrueacf(n = n, H = H)
        plot(ans, type = "h", col = "steelblue")
        title(main = "FARIMA True ACF")
        grid()
        abline(h = 0, col = "grey")

        # FGN FFT:
        ans = farimaTruefft(n = n, H = H)
        plot(Re(ans), type = "h", col = "steelblue")
        title(main = "FARIMA True FFT")
        grid()
        abline(h=0, col = "grey")

        # Reset Frame:
        par(mfrow = c(1, 1), cex = 0.7)
    }

    # Open Slider Menu:
    fBasics:::.sliderMenu(refresh.code,
       names =       c(  "n",    "H"),
       minima =      c(   10,   0.01),
       maxima =      c(  200,   0.99),
       resolutions = c(   10,   0.01),
       starts =      c(  100,   0.70))
}


# ------------------------------------------------------------------------------


.gkFARIMA0 =
function(n, H)
{
    # Description:
    #   Calculates  gk=fft of V=(r(0),...,r(n-2),r(n-1),r(n-2),...,r(1)),
    #   where r = the autocovariances of a fractional ARIMA with innovation
    #   variance 0

    # Arguments:
    #   n = length of time series
    #   H = self-similarity parameter

    # Value:
    #   gk = Fourier transform of V at the Fourier frequencies

    # Author:
    #   Jan Beran; modified: Martin Maechler, Date: Sep 95.

    # FUNCTION:

    # FFT:
    gammak = .ckFARIMA0(n,H)
    ind = c(0:(n - 2), (n - 1), (n - 2):1)
    gk = gammak[ind+1]
    ans = drop(fft(c(gk), inverse = TRUE))

    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


.simFARIMA0 =
function(n, H)
{
    # Description:
    #   Simulates a series X(1),...,X(n) of a fractional ARIMA(0,d,0)
    #   process (d=H-1/2)

    # Arguments:
    #   n = length of time series
    #   H = self-similarity parameter

    # Value:
    #   Simulated series X(1),...,X(n)

    # Author:
    #   Jan Beran; modified: Martin Maechler, Date: Sep 95.

    # FUNCTION:

    # Simulate:
    z = rnorm(2*n)
    zr = z[c(1:n)]
    zi = z[c((n+1):(2*n))]
    zic = -zi
    zi[1] = 0
    zr[1] = zr[1]*sqrt(2)
    zi[n] = 0
    zr[n] = zr[n]*sqrt(2)
    zr = c(zr[c(1:n)], zr[c((n-1):2)])
    zi = c(zi[c(1:n)], zic[c((n-1):2)])
    z = complex(real = zr, imaginary = zi)
    cat("n = ", n, "h = ", H)
    gksqrt = Re(.gkFARIMA0(n, H))
    if (all(gksqrt > 0)) {
        gksqrt = sqrt(gksqrt)
        z = z*gksqrt
        z = fft(z, inverse = TRUE)
        z = 0.5*(n-1)**(-0.5)*z
        ans = drop(Re(z[c(1:n)]))
    } else {
        gksqrt = 0*gksqrt
        stop("Re(gk)-vector not positive")
    }

    # Return Value:
    ans
}


################################################################################

