
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
# FUNCTIONS:            FRACTIONAL GAUSSIAN NOISE:
#  fgnSim                Generates fractional Gaussian noise
#  .fgnSim.durbin         Durbin's Method
#  .fgnSim.paxson         Paxson's Method
#  .fgnSim.beran          Beran's Method
#  .fgnSlider             Displays fgn simulated time Series
################################################################################


fgnSim =
function(n = 1000, H = 0.7, method = c("beran", "durbin", "paxson"))
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Creates a series of fractional Gaussian Noise

    # Arguments:
    #   n - length of desired series.
    #   H - the Hurst parameter.
    #   sigma - the standard deviation of the innovations used.

    # Details:
    #   FGN time series simulation. The FGN sequences use
    #   functions going back to Taqqu et al., Paxson and
    #   Beran (with Maechler's modifications from StatLib).

    # FUNCTION:

    # Settings:
    mean = 0
    std = 1
    method = match.arg(method)
    ans = NA

    # Generate Sequence:
    if (method == "beran")
        ans = .fgnSim.beran (n = n, H = H, mean = mean, std = std)
    if (method == "durbin")
        ans = .fgnSim.durbin(n = n, H = H, mean = mean, std = std)
    if (method == "paxson")
        ans = .fgnSim.paxson(n = n, H = H, mean = mean, std = std)
    if (is.na(ans[1])) stop("No valid method selected.")

    # Result:
    ans = as.ts(ans)
    attr(ans, "control") <- c(method = method, H = H)

    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


.fgnSim.durbin =
function(n, H, mean, std)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Function to simulate a FGN  sequence, using the Durbin-Levinson
    #   coefficients, along the original C Function of Vadim Teverovsky
    #   Original Cdurbin.C, double loop removed - now single loop in R!
    #   This makes the function quite fast.

    # Settings:
    n = n + 1
    h = H
    sigma = std
    normal = rnorm(n+1)
    sigma2 = sigma*sigma/2
    acov0 = 2*sigma2
    acov = vee = phi1 = phi2 = output = rep(0, n)
    I = 1:n
    acov = sigma2 *( (I+1)^(2*h) - 2*I^(2*h) + abs(I-1)^(2*h) )
    phi1[1] = acov[1]/acov0
    phi2[1] = phi1[1]
    vee0 = acov0
    vee[1] = vee0 * (1 - phi1[1]^2)
    output[1] = sqrt(vee0) * normal[1]

    # Durbin-Levinson:
    for (i in 2:n){
        phi1[i] = acov[i]
        J = 1:(i-1)
        phi1[i] = phi1[i] - sum(phi2[J]*acov[i-J])
        phi1[i] = phi1[i]/vee[i-1]
        vee[i] = vee[i-1]*(1-phi1[i]^2)
        output[i] = sqrt(vee[i-1]) * normal[i]
        phi1[J] = phi2[J] - phi1[i]*phi2[i-J]
        output[i] = output[i] + sum(phi2[J] * output[i-J])
        phi2[1:i] =  phi1[1:i]
    }

    # Result:
    ans = sigma*output[-1] + mean
    attr(ans, "control") <- c(method = "durbin", H = H)

    # Return value:
    ans
}


# ------------------------------------------------------------------------------


.fgnSim.paxson =
function(n, H, mean, std)
{
    # Description:
    #   Generates a FGN sequence by Paxson's FFT-Algorithm

    # Details:
    #   This file contains a function for synthesizing approximate
    #   fractional Gaussian noise.
    #   * Note that the mean and variance of the returned points is
    #     irrelevant, because any linear transformation applied to the
    #     points preserves their correlational structure (and hence
    #     their approximation to fractional Gaussian noise); and by
    #     applying a linear transformation you can transform the
    #     points to have any mean and variance you want.
    #   * If you're using the sample paths for simulating network
    #     arrival counts, you'll want them to all be non-negative.
    #     Hopefully you have some notion of the mean and variance
    #     of your desired traffic, and can apply the corresponding
    #     transformation. If, after transforming, a fair number of
    #     the points are still negative, then perhaps your traffic
    #     is not well-modeled using fractional Gaussian noise.
    #     You might instead consider using an exponential transformation.
    #   * All of this is discussed in the technical report:
    #     Fast Approximation of Self-Similar Network Traffic,
    #     Vern Paxson, technical report LBL-36750/UC-405, April 1995.
    #     URL:ftp://ftp.ee.lbl.gov/papers/fast-approx-selfsim.ps.Z

    # Value:
    #   FGN vector of length n.

    # FUNCTION:

    # Returns a Fourier-generated sample path of a "self similar" process,
    # consisting of n points and Hurst parameter H (n should be even).
    n = n/2
    lambda = ((1:n)*pi)/n

    # Approximate ideal power spectrum:
    d = -2*H - 1
    dprime = -2*H
    a = function(lambda,k) 2*k*pi+lambda
    b = function(lambda,k) 2*k*pi-lambda
    a1 = a(lambda,1); b1 = b(lambda,1)
    a2 = a(lambda,2); b2 = b(lambda,2)
    a3 = a(lambda,3); b3 = b(lambda,3)
    a4 = a(lambda,4); b4 = b(lambda,4)
    FGNBest = a1^d+b1^d+a2^d+b2^d+a3^d+b3^d +
        (a3^dprime+b3^dprime + a4^dprime+b4^ dprime)/(8*pi*H)
    f = 2 * sin(pi*H) * gamma(2*H+1) * (1-cos(lambda)) *
        (lambda^(-2*H-1) + FGNBest )

    # Adjust for estimating power:
    # spectrum via periodogram.
    f = f * rexp(n)

    # Construct corresponding complex numbers with random phase.
    z = complex(modulus = sqrt(f), argument = 2*pi*runif(n))

    # Last element should have zero phase:
    z[n] = abs(z[n])

    # Expand z to correspond to a Fourier:
    # transform of a real-valued signal.
    zprime = c(0, z, Conj(rev(z)[-1]))

    # Inverse FFT gives sample path:
    z = Re(fft(zprime, inverse = TRUE))

    # Standardize:
    z = (z-mean(z))/sqrt(var(z))
    ans = std*z + mean
    attr(ans, "control") <- c(method = "paxson", H = H)

    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


.fgnSim.beran =
function(n, H = 0.7, mean = 0, std = 1)
{
    # Description:
    #   Generates a FGN sequence by Beran's FFT-Algorithm

    # Value:
    #   FGN vector of length n.

    # FUNCTION:

    # Generate Sequence:
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
    z = complex(real = zr,imaginary = zi)

    # .gkFGN0:
    k = 0:(n-1)
    gammak = (abs(k-1)**(2*H)-2*abs(k)**(2*H)+abs(k+1)**(2*H))/2
    ind = c(0:(n - 2), (n - 1), (n - 2):1)
    .gkFGN0 = fft(c(gammak[ind+1]), inverse = TRUE)
    gksqrt = Re(.gkFGN0)
    if (all(gksqrt > 0)) {
        gksqrt = sqrt(gksqrt)
        z = z*gksqrt
        z = fft(z, inverse = TRUE)
        z = 0.5*(n-1)**(-0.5)*z
        z = Re(z[c(1:n)])
    } else {
        gksqrt = 0*gksqrt
        stop("Re(gk)-vector not positive")
    }

    # Standardize:
    # (z-mean(z))/sqrt(var(z))
    ans = std*drop(z) + mean
    attr(ans, "control") <- c(method = "beran", H = H)

    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


.fgnSlider =
function()
{   # A function implemented by Diethelm Wuertz

    # Description
    #   Displays fgn simulated time Series

    # FUNCTION:

    # Internal Function:
    refresh.code = function(...)
    {
        # Sliders:
        n      = fBasics:::.sliderMenu(no = 1)
        H      = fBasics:::.sliderMenu(no = 2)
        method = fBasics:::.sliderMenu(no = 3)

        # Graph Frame:
        par(mfrow = c(1, 1))

        # Select Method:
        Method = c("beran", "durbin", "paxson")
        Method = Method[method]

        # FGN TimeSeries:
        x = fgnSim(n = n, H = H, method = Method)
        plot(x, type = "l", col = "steelblue")
        grid()
        abline(h = 0, col = "grey")
        title(main = paste(Method, "FGN | H =", H))

        # Reset Frame:
        par(mfrow = c(1, 1))
    }

    # Open Slider Menu:
    fBasics:::.sliderMenu(refresh.code,
       names =       c(  "n",    "H", "method"),
       minima =      c(   10,   0.01,       1),
       maxima =      c(  200,   0.99,       3),
       resolutions = c(   10,   0.01,       1),
       starts =      c(  100,   0.70,       1))
}


################################################################################

