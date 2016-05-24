
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
#   1999 - 2008, Diethelm Wuertz, Rmetrics Foundation, GPL
#   Diethelm Wuertz <wuertz@itp.phys.ethz.ch>
#   info@rmetrics.org
#   www.rmetrics.org
# for the code accessed (or partly included) from other R-ports:
#   see R's copyright and license files
# for the code accessed (or partly included) from contributed R-ports
# and other sources
#   see Rmetrics's copyright file


################################################################################
# FUNCTION:               FORECASTING: 
#  garchKappa              Computes Expection for APARCH Models
#  .funE                   Internal function used by kappa()
################################################################################


test.predict.methods <-  
    function()
{   
    # RVs:
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")
    
    print(NA)
    
    # Return Value:
    return()    
} 
 

# ------------------------------------------------------------------------------


test.garchKappa <-  
    function()
{
    # garchKappa(
    #   cond.dist = c("norm", "ged", "std", "snorm", "sged", "sstd"), 
    #   gamma = 0, delta = 2, skew = NA, shape = NA)
    
    garchKappa()
    
    garchKappa("norm", gamma = 0.3)
    
    garchKappa("ged", skew = 0.95, shape = 1)
    
    # CHECK !!! garchKappa("dstd", skew = 0.95, shape = 1)
    
    # Return Value:
    return()    
} 
 

################################################################################
    
