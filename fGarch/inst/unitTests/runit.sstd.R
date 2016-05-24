
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
# FUNCTION:              VARIANCE-1 STUDENT-T DISTRIBUTION:
#  dstd                   Density for the Student-t Distribution
#  pstd                   Probability function for the Student-t Distribution
#  qstd                   Quantile function for the Student-t Distribution
#  rstd                   Random Number Generator for the Student-t
# FUNCTION:              SKEW VARIANCE-1 STUDENT-T DISTRIBUTION:
#  dsstd                  Density for the skewed Student-t Distribution
#  psstd                  Probability function for the skewed STD
#  qsstd                  Quantile function for the skewed STD
#  rsstd                  Random Number Generator for the skewed STD
#  stdSlider              Displays Variance-1 Student-t Distribution and RVS
# FUNCTION:              PARAMETER ESTIMATION:              
#  stdFit                 Fit the parameters for a Sudent-t distribution
#  sstdFit                Fit the parameters for a skew Sudent-t distribution
################################################################################


test.sstdDist <-  
    function()
{ 
    # Standardized Student-t Distribution:
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")
    
    # Test:
    test = .distCheck("std",  mean = 0, sd = 1, nu = 5, robust = FALSE) 
    print(test)                                    
    
    # Skew Standardized Student-t Distribution:
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")
    
    # Test:
    test = .distCheck("sstd", mean = 0, sd = 1, nu = 5, xi = 1.5, robust = FALSE) 
    print(test)                                     
    
    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.stdFit <- 
    function()
{         
    # Fit the parameters for a Student-t distribution
    # stdFit - Fit the parameters for a Sudent-t distribution
    
    # Standardized Student-t Distribution:
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")
    
    # Series:
    x = rstd(n = 2500, mean = 0, sd = 1, nu = 5)
    
    # Fit:
    fit = stdFit(x)
    print(fit)
             
    # Fit the parameters for a skew Sudent-t distribution
    # sstdFit - Fit the parameters for a Sudent-t distribution
    
    # Skew Standardized Student-t Distribution:
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")
    
    # Series:
    x = rsstd(n = 2500, mean = 0, sd = 1, nu = 5, xi = 1.5)
    
    # Fit:
    fit = sstdFit(x)
    print(fit)
    
    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.sstdSlider <- 
    function()
{   
    # Try Distribution:
    # sstdSlider(type = "dist")                              
    NA
    
    # Try Random Variates:
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")
    
    # sstdSlider(type = "rand")
    NA
    
    # Return Value:
    return()    
}



################################################################################
    
