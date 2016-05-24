
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
# FUNCTION:           FRACTIONAL BROWNIAN MOTION:
#  fbmSim              Generates fractional Brownian motion
#  fgnSim              Generates fractional Gaussian noise
#  farimaSim           Generates FARIMA time series process
# FUNCTION :          FROM BERAN'S SPLUS SCRIPTS:
#  whittleFit          Whittle estimator
# FUNCTION:           HURST EXPONENT A LA TAQQU ET AL:
#  fHURST              S4 Class Representation
#   print.fHURST        S3 Print Method
#   plot.fHURST         S3 Plot Method
#  aggvarFit           3.1 Aggregated variance method
#  diffvarFit          3.2 Differenced aggregated variance method
#  absvalFit           3.3 Absolute values (moments) method
#  higuchiFit          3.4 Higuchi's method
#  pengFit             3.5 peng's or Residuals of Regression method
#  rsFit               3.6 R/S method
#  perFit              3.7 Periodogram and cumulated periodogram method
#  boxperFit           3.8 Boxed (modified) peridogram method
#  whittleFit          3.9 Whittle estimator -> PART II
#  hurstSlider         Hurst Slider [1-7]
# FUNCTIONS:          WAVELET ESTIMATOR:
#  waveletFit          Wavelet Estimator
################################################################################


test.fbmSim =
function()
{
    # Simulation:

    # Try:
    #   .fbmSimSlider()

    # fbmSim(
    #   n = 100, H = 0.7,
    #   method = c("mvn", "chol", "lev", "circ", "wave"),
    #   waveJ = 7,
    #   doplot = TRUE, fgn = FALSE)

    # Graphics Frame:
    par(mfrow = c(3, 2), cex = 0.7)

    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")
    x = fbmSim(n = 50, method = "mvn")
    print(x)
    target = round(mean(x), 2)
    print(target)
    current = +0.05
    checkEquals(target, current)

    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")
    x = fbmSim(n = 50, method = "chol")
    print(x)
    target = round(mean(x), 2)
    print(target)
    current = +0.45
    checkEquals(target, current)

    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")
    x = fbmSim(n = 50, method = "lev")
    print(x)
    target = round(mean(x), 2)
    print(target)
    current = 0.66
    checkEquals(target, current)

    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")
    x = fbmSim(n = 50, method = "circ")
    print(x)
    target = round(mean(x), 2)
    print(target)
    current = -0.45
    checkEquals(target, current)

    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")
    x = fbmSim(n = 50, method = "wave")
    print(x)
    target = round(mean(x), 2)
    print(target)
    current = 0.09
    checkEquals(target, current)

    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.fgnSim =
function()
{
    # Simulation:

    # Try:
    #   .fgnSimSlider()

    # fbmSim(
    #   n = 100, H = 0.7,
    #   method = c("mvn", "chol", "lev", "circ", "wave"),
    #   waveJ = 7, doplot = TRUE, fgn = FALSE)

    # Graphics Frame:
    par(mfrow = c(3, 2), cex = 0.7)

    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")


    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.farimaSim =
function()
{
    # Simulation:

    # Try:
    #   .farimaSimSlider()

    # farimaSim(
    #   n = 1000,
    #   model = list(ar = c(0.5, -0.5), d = 0.3, ma = 0.1),
    #   method = c("freq", "time"),
    #   ...)

    # Frequency Method:
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")
    x = farimaSim(n = 50, method = "freq")
    print(x)
    target = round(var(x), 2)
    print(target)
    current = 1.6
    print(current)
    checkEquals(target, current)

    # Time Methhod:
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")
    x = farimaSim(n = 50, method = "time")
    print(x)
    target = round(var(x), 2)
    print(target)
    current = 1.6
    print(current)
    checkEquals(target, current)

    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


## test.whittleFit =
## function()
## {
##     # whittleFit(
##     #   x, order = c(1, 1), subseries = 1,
##     #   method = c("fgn", "farma"),
##     #   trace = FALSE, spec = FALSE, title = NULL, description = NULL)


##     # Beran - Simulate:
##     RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
##     set.seed(4711, kind = "Marsaglia-Multicarry")
##     x = fgnSim(n = 1000, H = 0.7)

##     # Fit:
##     Hurst = whittleFit(x)@hurst$H
##     print(Hurst)
##     target = round(Hurst, 3)
##     print(target)
##     current = 0.701
##     checkEquals(target, current)

##     # Durbin - Simulate:
##     RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
##     set.seed(4711, kind = "Marsaglia-Multicarry")
##     x = fgnSim(n = 1000, H = 0.7, method = "durbin")

##     # Fit:
##     Hurst = whittleFit(x)@hurst$H
##     print(Hurst)
##     target = round(Hurst, 3)
##     print(target)
##     current = 0.679
##     checkEquals(target, current)

##     # Paxson - Simulate:
##     RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
##     set.seed(4711, kind = "Marsaglia-Multicarry")
##     x = fgnSim(n = 1000, H = 0.7, method = "paxson")

##     # Fit:
##     Hurst = whittleFit(x)@hurst$H
##     print(Hurst)
##     target = round(Hurst, 3)
##     print(target)
##     current = 0.728
##     checkEquals(target, current)

##     # Return Value:
##     return()
## }


# ------------------------------------------------------------------------------


test.hurstFit =
function()
{
    # Try:
    #   hurstSlider(x = fbmSim(n = 1000, H = 0.7, fgn = TRUE))

    # Methods:
    #   method = 1: aggvarFit
    #   method = 2: diffvarFit
    #   method = 3: absvalFit
    #   method = 4: higuchiFit
    #   method = 5: pengFit
    #   method = 6: rsFit
    #   method = 7: perFit

    # Simulate:
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")
    x = fgnSim(n = 1000, H = 0.7)

    Hurst = aggvarFit(x)@hurst$H
    print(Hurst)
    target = round(Hurst, 3)
    print(target)
    current = 0.731
    checkEqualsNumeric(target, current)

    Hurst = diffvarFit(x)@hurst$H
    print(Hurst)
    target = round(Hurst, 3)
    print(target)
    current = 0.869
    checkEqualsNumeric(target, current)

    Hurst = absvalFit(x)@hurst$H
    print(Hurst)
    target = round(Hurst, 3)
    print(target)
    current = 0.787
    checkEqualsNumeric(target, current)

    Hurst = higuchiFit(x)@hurst$H
    print(Hurst)
    target = round(Hurst, 3)
    print(target)
    current = 0.714
    checkEqualsNumeric(target, current)

    Hurst = pengFit(x)@hurst$H
    print(Hurst)
    target = round(Hurst, 3)
    print(target)
    current = 0.620
    checkEqualsNumeric(target, current)

    Hurst = rsFit(x)@hurst$H
    print(Hurst)
    target = round(Hurst, 3)
    print(target)
    current = 0.594
    checkEqualsNumeric(target, current)

    Hurst = perFit(x)@hurst$H
    print(Hurst)
    target = round(Hurst, 3)
    print(target)
    current = 0.602
    checkEqualsNumeric(target, current)

    # More Estimators:
    # ...

    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.waveletFit =
function()
{
    # waveletFit(
    #   x, length = NULL, order = 2, octave = c(2, 8),
    #   doplot = FALSE, title = NULL, description = NULL)

    # Simulate:
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")
    x = fgnSim(n = 1000, H = 0.7)

    # Fit:
    Hurst = waveletFit(x)@hurst$H
    print(Hurst)
    target = round(Hurst, 3)
    print(target)
    current = 0.696
    checkEquals(target, current)

    # Return Value:
    return()
}


################################################################################
