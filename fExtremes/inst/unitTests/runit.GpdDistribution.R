
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
# FUNCTION:              GPD DISTRIBUTION FAMILY:
# dgpd                    Density for the Generalized Pareto DF [USE FROM EVIS]
#  pgpd                    Probability for the Generalized Pareto DF
#  qgpd                    Quantiles for the Generalized Pareto DF
#  rgpd                    Random variates for the Generalized Pareto DF
# gpdMoments              Computes true statistics for GPD distribution
# gpdSlider               Displays distribution and rvs for GPD distribution
################################################################################


test.gpd =
function()
{
    # Check Distribution:
    set.seed(1985)
    .distCheck(fun = "gpd", n = 500, xi = 1, mu = 0, beta = 1)
    
    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.gpdMoments = 
function()
{
    # gpdMoments(xi = 1, mu = 0, beta = 1) 

    # Compute Moments:
    xi = seq(-2, 2, length = 401)
    mom = gpdMoments(xi)
    
    # Plot Mean:
    par(mfrow = c(2, 1), cex = 0.7)
    par(ask = FALSE)
    plot(xi, mom$mean, main = "Mean", pch = 19, cex = 0.5)
    abline(v = 1, col = "red", lty = 3)
    abline(h = 0, col = "red", lty = 3)
    
    # Plot Variance:
    plot(xi, log(mom$var), main = "log Variance", pch = 19, cex = 0.5)
    abline(v = 1/2, col = "red", lty = 3)
    abline(h = 0.0, col = "red", lty = 3)
    
    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.gpdSlider = 
function()
{
    # Distribution Slider:
    # print("Activate Slider manually!")
    # gpdSlider(method = "dist")
    
    # Random Variates Slider:
    # print("Activate Slider manually!")
    # gpdSlider(method = "rvs")
    NA
    
    # Return Value:
    return()    
}


################################################################################

