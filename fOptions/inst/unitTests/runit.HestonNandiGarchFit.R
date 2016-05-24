
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
#  hngarchSim            Simulates an HN-GARCH(1,1) Time Series Process
#  hngarchFit            Fits a HN-GARCH model by Gaussian Maximum Likelihood
#  print.hngarch         Print method, reports results
#  summary.hngarch       Summary method, diagnostic analysis
#  hngarchStats          Computes Unconditional Moments of a HN-GARCH Process
################################################################################


test.hngarchSim = 
function()
{
    # Simulate a Heston-Nandi Garch(1,1) Process
    
    # RVs:
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")
    
    # Symmetric Model - Parameters:
    model = list(lambda = 4, omega = 8e-5, alpha = 6e-5, 
        beta = 0.7, gamma = 0, rf = 0)
        
    # Series:
    x = hngarchSim(model = model, n = 500, n.start = 100)
    
    # Plot:
    par(mfrow = c(2, 1), cex = 0.75)
    plot(x, type = "l", col = "steelblue", main = "HN Garch Symmetric Model")
    grid()

    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.hngarchFit = 
function()
{    
    # Simulate a Heston-Nandi Garch(1,1) Process:
    
    # RVs:
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")
    
    # Symmetric Model - Parameters:
    model = list(lambda = 4, omega = 8e-5, alpha = 6e-5, 
        beta = 0.7, gamma = 0, rf = 0)
    x = hngarchSim(model = model, n = 500, n.start = 100)
                                                                     
    # Estimate Parameters:
    # HN-GARCH log likelihood Parameter Estimation:
    # To speed up, we start with the simulated model ...
    
    # Fit Symmetric Case:
    mle = hngarchFit(x = x, model = model, trace = TRUE, symmetric = TRUE)
    print(mle)
    
    # Assymmetric Case:
    mle = hngarchFit(x = x, model = model, trace = TRUE, symmetric = FALSE)
    print(mle)
        
    # HN GARCH Plot:
    # ... there is no plot - plotting is done in summary 
    
    # HN-GARCH Diagnostic Analysis:
    # Note, residuals are still missing ...
    par(mfrow = c(3, 1))
    summary(mle, col = "steelblue")                                             
    
    # HN-GARCH Moments:
    hngarchStats(mle$model)    

    # Return Value:
    return()    
}


################################################################################

