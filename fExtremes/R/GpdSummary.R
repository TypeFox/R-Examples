
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
# METHODS:                PRINT, PLOT, AND SUMMARY:
#  summary.fGPDFIT         S3 Summary Method for object of class "fGPDFIT"
################################################################################


summary.fGPDFIT = 
function(object, doplot = TRUE, which = "all", ...) 
{   # A function written by Diethelm Wuertz

    # Description:
    #   Summary method for objects of class "gpdFit"
    
    # Arguments:
    
    # FUNCTION:

    # Title:
    cat("\nTitle:\n" , object@title, "\n")
    
    # Function Call:
    cat("\nCall:\n")
    cat(paste(deparse(object@call), sep = "\n", 
        collapse = "\n"), "\n", sep = "") 
            
    # Estimation Type:
    cat("\nEstimation Type:\n ", object@method, "\n") 
    
    # Estimated Parameters:
    cat("\nEstimated Parameters:\n")
    print(object@fit$par.ests)
    
    # Summary - For MLE print additionally:
    if (object@method[2] == "mle") {
        cat("\nStandard Deviations:\n"); print(object@fit$par.ses)
        if (!is.na(object@fit$llh))
            cat("\nLog-Likelihood Value:\n ", object@fit$llh, "\n")
        if (!is.na(object@fit$convergence))
            cat("\nType of Convergence:\n ", object@fit$convergence, "\n") 
    }
    
    # Plot:
    if (doplot) {
        plot(object, which = which, ...)
    }
    
    # Desription:
    cat("\nDescription\n ", object@description, "\n\n")
    
    # Return Value:
    invisible(object)
}


################################################################################

