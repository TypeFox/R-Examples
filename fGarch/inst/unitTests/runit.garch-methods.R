
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
# METHODS:                DESCRIPTION:
#  show.fGARCH             S4 print method for an object of class 'fGARCH'
#  summary.fGARCH          S3 summary method for an object of class 'fGARCH'
#  plot.fGARCH             S3 plot method for an object of class 'fGARCH'
#  residuals.fGARCH        S3 residuals method for an object of class 'fGARCH'
#  fitted.fGARCH           S3 fitted values for an object of class 'fGARCH'
#  predict.fGARCH          S3 prediction method for an object of class 'fGARCH'
# STATISTICS:             Description:
#  .truePersistence        Compute persistence
################################################################################


test.garch.methods.show <- 
    function()
{ 
    # show.fGARCH - S4 print method for an object of class 'fGARCH'
    
    # Garch(1,1) Default Model:
    x = garchSim(n = 250)
    
    # Fit:
    fit = garchFit( ~ garch(1,1), data = x, trace = FALSE)
    
    # Note:
    # Does also work:
    #   fit = garchFit(dem2gbp ~garch(1,1), data = dem2gbp, trace = FALSE)
    # Note:
    #   fit = garchFit(DEM2GBP ~garch(1,1), data = dem2gbp, trace = FALSE)
    #   [1] "DEM2GBP"
    #   [1] "dem2gbp"
    #   Error in .garchArgsParser(formula = formula, data = data, trace = FALSE) : 
    #     Formula and data units do not match.

    # Print:
    print(fit)
    show(fit)
    
    # Summary:
    summary(fit)
    
    # Return Value:
    return()    
} 


################################################################################
    
