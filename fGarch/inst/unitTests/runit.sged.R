
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
# FUNCTION:              GED DISTRIBUTION:
#  dged                   Density for the Generalized Error Distribution
#  pged                   Probability function for the GED
#  qged                   Quantile function for the GED
#  rged                   Random Number Generator for the GED
# FUNCTION:              SKEW GED DISTRIBUTION:
#  dsged                  Density for the skewed GED
#  psged                  Probability function for the skewed GED
#  qsged                  Quantile function for the skewed GED
#  rsged                  Random Number Generator for the skewed GED
# FUNCTION:              GED DISTRIBUTION SLIDER:
#  sgedSlider            Displays Generalized Error Distribution and RVS
# FUNCTION:              PARAMETER ESTIMATION:
#  gedFit                 Fit the parameters for a GED distribution
#  sgedFit                Fit the parameters for a skew GED distribution
# FUNCTION:              MOMENTS:
################################################################################


test.sgedDis <- 
    function()
{       
    # Generalized Error Distribution:
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")
    
    # Test:
    test = .distCheck("ged",  mean = 0, sd = 1, nu = 2, robust = FALSE) 
    print(test)                                     
       
    # Skew Generalized Error Distribution:
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(1953, kind = "Marsaglia-Multicarry")
    test = .distCheck("sged", mean = 0, sd = 1, nu = 2, xi = 0.8, robust = FALSE) 
    print(test)                                       
    
    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.sgedFit <- 
    function()
{      
    # Parameter Estimation:
    #  gedFit - Fit the parameters for a GED distribution
    
    # Generalized Error Distribution:
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")
    
    # Series:
    x = rged(1000, mean = 0, nu = 2)
    
    # Fit:
    fit = gedFit(x)
    print(fit)
      
    # Fit the parameters for a skew GED distribution
    #  sgedFit - Fit the parameters for a skew GED distribution
    
    # Skew Generalized Error Distribution:
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")
    x = rsged(1000, mean = 0, sd = 1, nu = 2, xi = 1.5)
    fit = sgedFit(x)
    print(fit)
    
    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.sgedSlider <- 
    function()
{   
    # Try Distribution:
    # sgedSlider(type = "dist")
    NA 
    
    # Try Random Variates:
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")
    
    # sgedSlider(type = "rand")
    NA
    
    # Return Value:
    return()    
}


################################################################################
    
