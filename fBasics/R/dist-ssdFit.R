
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


################################################################################
# FUNCTIONS:       DESCRIPTION
#   ssdFit          Estimate probability densities using smoothing spline ANOVA
#   .print.ssd      S3 Print Method
################################################################################


ssdFit <-  
    function(x)
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #    Estimate probability densities using smoothing spline ANOVA 
    #    models with cubic spline, linear spline, or thin-plate spline 
    #    marginals for numerical variables.
    
    # FUNCTION:
    
    # Parameter Fit: 
    ans <- gss::ssden(~x)
    
    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


.print.ssd <- 
    function(x, ...)
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    
    # FUNCTION:
    
    # call
    cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
    
    # Terms:
    cat("Terms:\n")
    print.default(x$terms$labels)
    cat("\n")
    
    # terms overview
    cat("Number of unpenalized and penalized terms:\n\n")
    print.default(x$desc)
    cat("\n")
    cat("Smoothing parameters are selected by CV with alpha=", x$alpha, ".", 
        sep = "")
    cat("\n")
    
    # the rest are suppressed
    invisible()
}


################################################################################

