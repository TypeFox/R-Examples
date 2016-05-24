
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR Description. See the 
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General 
# Public License along with this library; if not, write to the 
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston, 
# MA 02111-1307 USA


################################################################################
# FUNCTION:                     DESCRIPTION:             
#  summary.fPORTFOLIO            S3 Summary method for 'fPORTFOLIO' objects
################################################################################


summary.fPORTFOLIO <- 
    function(object, ...)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Plot method for an object of class 'fPORTFOLIO'
    
    # Note:
    #   This method can also be used for plotting graphs fitted by 
    #   the function 'garch' from the contributed R package 'tseries'.
    
    # FUNCTION:

    # Summary:
    print(object)
    funCalled = as.character(object@call[1])
    if (funCalled == "portfolioFrontier") {      
        weightsPlot(object)
        weightedReturnsPlot(object)
        covRiskBudgetsPlot(object)
        # Plot Frontier:
        plot(object, which = 1)
    } else {
        weightsPie(object)
        weightedReturnsPie(object)
        covRiskBudgetsPie(object)
    }
          
    # Return Value:
    invisible(object)
} 


################################################################################

