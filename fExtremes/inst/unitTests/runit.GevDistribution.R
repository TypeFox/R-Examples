
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
# FUNCTION:             GEV DISTRIBUTION FAMILY: [CALLING EVD]
#  dgev                  Density for the GEV Distribution
#   pgev                  Probability for the GEV Distribution
#   qgev                  Quantiles for the GEV Distribution
#   rgev                  Random variates for the GEV Distribution
#  gevMoments            Computes true statistics for GEV distribution
#  gevSlider             Displays distribution and rvs for GEV distribution
################################################################################


test.gev =
function()
{
    # Check Distribution:
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")
    .distCheck(fun = "gev", n = 2000, xi = 0.0, mu = 0, beta = 1)

    # Check Distribution:
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")
    .distCheck(fun = "gev", n = 5000, xi = 0.3, mu = 0, beta = 2)

    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.gevMoments =
function()
{
    # gevMoments(xi = 0, mu = 0, beta = 1)

    # Compute Moments:
    xi = seq(-4.5, 1.5, by = 0.25)
    mom = gevMoments(xi)
    print(mom)

    # Plot Mean:
    par(mfrow = c(2, 1), cex = 0.7)
    par(ask = FALSE)
    xi = seq(-5, 2, length = 351)
    mom = gevMoments(xi)
    plot(xi, mom$mean, main = "Mean GEV", pch = 19, col = "steelblue")
    abline(v = 1, col = "red", lty = 3)
    abline(h = 0, col = "red", lty = 3)

    # Plot Variance:
    plot(xi, log(mom$var), main = "log Variance GEV", pch = 19, col = "steelblue")
    abline(v = 1/2, col = "red", lty = 3)
    abline(h = 0.0, col = "red", lty = 3)


    # check gevMoments for specific values
    xi <- c(-1, 0, 0.3)
    mu <- c(-1, 0, 1)
    beta <- c(0.5, 1, 10)

    for (i in seq(length(xi))) {
        for (j in seq(length(xi))) {
            for (k in seq(length(xi))) {
                rg <- rgev(1000000, xi = xi[i], mu = mu[j], beta = beta[k])
                rgMoments <- gevMoments(xi = xi[i], mu = mu[j], beta = beta[k])
                checkEqualsNumeric(mean(rg), rgMoments$mean, tolerance = 0.1)
                checkEqualsNumeric(var(rg), rgMoments$var, tolerance = 0.1)
            }
        }
    }

    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.gevSlider =
function()
{
    # Distribution Slider:
    # print("Activate Slider manually!")
    # gevSlider(method = "dist")
    NA

    # Random Variates Slider:
    # print("Activate Slider manually!")
    # gevSlider(method = "rvs")
    NA

    # Return Value:
    return()
}


################################################################################

