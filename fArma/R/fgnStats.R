
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
#  fgnTrueacf            Returns FGN true autocorrelation function
#  fgnTruefft            Returns FGN true fast Fourier transform
#  .fgnStatsSlider        Displays fgn true statistics
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


fgnTrueacf =
function(n =100, H = 0.7)
{   # A function implemented by Diethelm Wuertz

    # FUNCTION:

    # True ACF:
    ans = .ckFGN0(n = n, H = H)

    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


fgnTruefft =
function(n = 100, H = 0.7)
{   # A function implemented by Diethelm Wuertz

    # FUNCTION:

    # True FFT:
    ans = .gkFGN0(n = n, H = H)

    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


fgnStatsSlider =
function()
{   # A function implemented by Diethelm Wuertz

    # Description
    #   Displays FGN true statistics: ACF and FFT

    # Example:
    #   fgnStatsSlider()

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
        ans = fgnTrueacf(n = n, H = H)
        plot(ans, type = "h", col = "steelblue")
        title(main = "FGN True ACF")
        grid()
        abline(h = 0, col = "grey")

        # FGN FFT:
        ans = fgnTruefft(n = n, H = H)
        plot(Re(ans), type = "h", col = "steelblue")
        title(main = "FGN True FFT")
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


.ckFGN0 =
function(n, H)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Computes covariances of a fractional Gaussian  process

    # Arguments:
    #   n = length of time series
    #   H = self-similarity parameter

    # Value:
    #   Covariances upto lag n-1

    # Author:
    #   Jan Beran; modified: Martin Maechler, Date: Sep 95.

    # FUNCTION:
    k = 0:(n-1)
    H2 = 2 * H
    ans = drop((abs(k-1)**H2-2*abs(k)**H2+abs(k+1)**H2)/2)

    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


.gkFGN0 =
function(n, H)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Calculates gk = fft of V=(r(0),...,r(n-2),
    #   r(n-1), r(n-2),...,r(1), where r=the autocovariances
    #   of a fractional Gaussian process with variance 1

    # Arguments:
    #   n = length of time series
    #   H = self-similarity parameter

    # Value:
    #   gk = Fourier transform of V at Fourier frequencies

    # Author:
    #   Jan Beran; modified: Martin Maechler, Date: Sep 95.

    # FUNCTION:

    # FFT:
    gammak = .ckFGN0(n, H)
    ind = c(0:(n - 2), (n - 1), (n - 2):1)
    gk = gammak[ind+1]
    ans = drop(fft(c(gk), inverse = TRUE))

    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


.simFGN0 =
function(n, H)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Simulates a series X(1),...,X(n) of a fractional Gaussian process

    # Arguments:
    #   n = length of time series
    #   H = self-similarity parameter

    # Value:
    #   simulated series X(1),...,X(n)

    # Author:
    #   Jan Beran; modified: Martin Maechler, Date: Sep 95.

    # FUNCTION:

    # Simulation:
    z = rnorm(2*n)
    zr = z[c(1:n)]
    zi = z[c((n+1):(2*n))]
    zic = -zi
    zi[1] = 0
    zr[1] = zr[1]*sqrt(2)
    zi[n] = 0
    zr[n] = zr[n]*sqrt(2)
    zr = c(zr[c(1:n)],zr[c((n-1):2)])
    zi = c(zi[c(1:n)],zic[c((n-1):2)])
    z = complex(real = zr,imaginary = zi)
    gksqrt = Re(.gkFGN0(n, H))
    if (all(gksqrt > 0)) {
        gksqrt = sqrt(gksqrt)
        z = z*gksqrt
        z = fft(z, inverse = TRUE)
        z = 0.5*(n-1)**(-0.5)*z
        ans = drop(Re(z[c(1:n)]))
    } else {
        gksqrt = 0*gksqrt
        cat("Re(gk)-vector not positive") }

    # Return Value:
    ans
}


################################################################################

