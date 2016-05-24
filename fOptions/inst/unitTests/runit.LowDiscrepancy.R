
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
# FUNCTION:             DESCRIPTION:
#  runif.pseudo           Uniform Pseudo Random number sequence
#  rnorm.pseudo           Normal Pseudo Random number sequence
#  runif.halton           Uniform Halton low discrepancy sequence
#  rnorm.halton           Normal Halton low discrepancy sequence
#  runif.sobol            Uniform Sobol low discrepancy sequence
#  rnorm.sobol            Normal Sobol low discrepancy sequence
################################################################################


test.pseudo =
function()
{
    # Pseudo Random Numbers:
    #   Uniform and Normal pseudo random number sequences

    # RVs:
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")

    # Graphics Frame:
    par(mfrow = c(2, 2), cex = 0.75)

    # Histogram Uniform:
    runif.pseudo(n = 10, dimension = 5)
    r = runif.pseudo(n = 1000, dimension = 1)
    hist(r, probability = TRUE, main = "Uniform Pseudo", xlab = "x",
        col = "steelblue", border = "white")
    abline (h = 1, col = "orange", lwd = 2)

    # Scatterplot Uniform:
    r = runif.pseudo(n = 1000, dimension = 2)
    plot(r, cex = 0.5, main = "Scatterplot Uniform Pseudo")

    # Histogram Normal:
    rnorm.pseudo(n = 10, dimension = 5)
    r = rnorm.pseudo(n = 1000, dimension = 1)
    hist(r, probability = TRUE, xlim = c(-3, 3), main = "Normal Pseudo",
        xlab = "x", col = "steelblue", border = "white")
    x = seq(-3, 3, length = 301)
    lines(x, dnorm(x), col = "orange", lwd = 2)

    # Scatterplot Normal:
    r = rnorm.pseudo(n = 1000, dimension = 2)
    plot(r, cex = 0.5, main = "Scatterplot Normal Pseudo")

    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.halton =
function()
{
    # Halton Sequence:
    #   Uniform and Normal Halton low discrepancy sequences

    # RVs:
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")

    # Graphics Frame:
    par(mfrow = c(2, 2), cex = 0.75)

    # Histogram Uniform:
    runif.halton(n = 10, dimension = 5)
    r = runif.halton(n = 5000, dimension = 1)
    hist(r, probability = TRUE, main = "Uniform Halton", xlab = "x",
        col = "steelblue", border = "white")
    abline (h = 1, col = "orange", lwd = 2)

    # Scatterplot Uniform:
    r = runif.halton(n = 1000, dimension = 2)
    plot(r, cex = 0.5, main = "Scatterplot Uniform Halton")

    # Histogram Normal:
    rnorm.halton(n = 10, dimension = 5)
    r = rnorm.halton(n = 5000, dimension = 1)
    hist(r, probability = TRUE, xlim = c(-3, 3), main = "Normal Halton",
        xlab = "x", col = "steelblue", border = "white")
    x = seq(-3, 3, length = 301)
    lines(x, dnorm(x), col = "orange", lwd = 2)

    # Scatterplot Normal:
    r = rnorm.halton(n = 1000, dimension = 2)
    plot(r, cex = 0.5, main = "Scatterplot Normal Halton")

    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.sobol =
function()
{
    # Sobol Sequence:
    #   Uniform and Normal Sobol low discrepancy sequences

    # RVs:
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")

    # Graphics Frame:
    par(mfrow = c(2, 2), cex = 0.75)

    # Histogram Uniform:
    runif.sobol(n = 10, dimension = 5)
    r = runif.sobol(5000, 1)
    hist(r, probability = TRUE, main = "Uniform Sobol",
        xlab = "x", col = "steelblue", border = "white")
    abline (h = 1, col = "orange", lwd = 2)

    # Scatterplot Uniform:
    r = runif.sobol(n = 1000, dimension = 2)
    plot(r, cex = 0.5, main = "Scatterplot Uniform Sobol")

    # Histogram Normal:
    rnorm.sobol(n = 10, dimension = 5)
    r = rnorm.sobol(1000, 1)
    hist(r, probability = TRUE, main = "Normal Sobol",
        xlab = "x", col = "steelblue", border = "white")
    x = seq(-3, 3, length = 301)
    lines(x, dnorm(x), col = "orange", lwd = 2)

    # Scatterplot Normal:
    r = rnorm.sobol(n = 1000, dimension = 2)
    plot(r, cex = 0.5, main = "Scatterplot Normal Sobol")
    
    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.scrambling =
function()
{
    # Sobol Scrambling:

    # RVs:
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")

    # runif.sobol(n, dimension, init = TRUE, scrambling = 0, seed = 4711)

    # Unscrambled:
    runif.sobol(10, 5)

    # Owen Type Scrambling:
    runif.sobol(10, 5, scrambling = 1)

    # Faure-Tezuka  Type Scrambling:
    runif.sobol(10, 5, scrambling = 2)

    # Combined Owen and Faure-Tezuka Type Scrambling:
    runif.sobol(10, 5, scrambling = 3)

    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.restart =
function()
{
    # Sobol Restart:

    # RVs:
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")

    # runif.sobol(n, dimension, init = TRUE, scrambling = 0, seed = 4711)
    runif.sobol(10, 5, init = TRUE)
    runif.sobol(10, 5, init = FALSE)

    # Seed:
    print(.getfOptionsEnv(".runif.sobol.seed"))

    # Return Value:
    return()
}


################################################################################

