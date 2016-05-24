
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
# FUNCTION:               GPD SIMULATION:
#  gpdSim                  Simulates a GPD distributed process
################################################################################


gpdSim = 
function(model = list(xi = 0.25, mu = 0, beta = 1), n = 1000, seed = NULL)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Generates random variates from a GPD distribution
    
    # Arguments:
    #   model - a list of model parameters, xi, mu and beta.
    #   n - an integer, the number of simulated random variates
    #   seed - an integer, the random number generator seed
    
    # FUNCTION:
    
    # Seed:
    if (is.null(seed)) seed = NA else set.seed(seed)
    
    # Simulate:
    ans = rgpd(n = n, xi = model$xi, mu = model$mu, beta = model$beta)  
    # DW: ans = as.ts(ans)
    ans = timeSeries(ans, units = "GPD")

    # Control:
    attr(ans, "control") = 
        data.frame(t(unlist(model)), seed = seed, row.names = "")
        
    # Return Value:
    ans 
}


################################################################################

