
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
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
# FUNCTION:                 DESCRIPTION:
#  show.fGARCH               S4 Show method for an object of class 'fGARCH'
#  show.fGARCHSPEC           S4 Show method for an object of class 'fGARCHSPEC'
################################################################################


setMethod(f = "show", signature(object = "fGARCH"), definition =
    function(object)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Print method for an object of class "fGARCH"

    # Arguments:
    #   object - an object of class 'fGARCH'

    # FUNCTION:

    # Title:
    cat("\nTitle:\n ")
    cat(object@title, "\n")

    # Call:
    cat("\nCall:\n ")
    cat(paste(deparse(object@call), sep = "\n", collapse = "\n"), "\n")

    # Mean and Variance Equation:
    cat("\nMean and Variance Equation:\n ")
    Name = unclass(attr(object@formula, "data"))
    Formula = object@formula
    attr(Formula, "data") <- NULL
    print(Formula)
    cat(" [", Name, "]\n", sep = "")

    # Univariate or Multivariate Modeling ?
    if (as.character(object@call[1]) == ".gogarchFit") 
    {
        # For multivariate Garch Models ...
        #   extract information from first fitted instrument.
        object@fit[[1]]@fit$params$cond.dist
        
        # Conditional Distribution:
        cat("\nConditional Distribution:\n ")
        cat(object@fit[[1]]@fit$params$cond.dist, "\n")
        
        # Number of Margins:
        cat("\nNumber of Margins:\n ")
        cat(length(object@fit), "\n")
        
    } else {
    
        # For univariate Garch Models ...
        
        # Conditional Distribution:
        cat("\nConditional Distribution:\n ")
        cat(object@fit$params$cond.dist, "\n")
    
        # Coefficients:
        cat("\nCoefficient(s):\n")
        digits = max(5, getOption("digits") - 4)
        print.default(format(object@fit$par, digits = digits), print.gap = 2,
             quote = FALSE)
    
        # Std. Errors:
        cat("\nStd. Errors:\n ")
        if (object@fit$params$cond.dist == "QMLE")
        {
            cat("robust", "\n")
        } else {
            cat("based on Hessian", "\n")
        }
        
        # Error Analysis:
        digits = max(4, getOption("digits") - 5)
        fit = object@fit
        signif.stars = getOption("show.signif.stars")
        cat("\nError Analysis:\n")
        printCoefmat(fit$matcoef, digits = digits, signif.stars = signif.stars)
    
        # Log Likelihood:
        cat("\nLog Likelihood:\n ")
        LLH = - object@fit$value
        N = NROW(object@data)
        cat(LLH, "   normalized: ", LLH/N, "\n")

    }
    
    # Description:
    cat("\nDescription:\n ")
    cat(object@description, "\n")

    # Return Value:
    cat("\n")
    invisible()
})


# ------------------------------------------------------------------------------


setMethod(f = "show", signature(object = "fGARCHSPEC"), definition =
    function(object)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   S4 Print Method for objects of class 'fGARCHSPEC'

    # Arguments:
    #   object - Object of class 'fGARCHSPEC'

    # FUNCTION:

    # Formula:
    x = object
    cat("\nFormula: \n ")
    cat(as.character(x@formula))

    # Model:
    cat("\nModel:")
    if (sum(abs(x@model$ar)) != 0)
        cat("\n ar:   ", x@model$ar)
    if (sum(abs(x@model$ma)) != 0)
        cat("\n ma:   ", x@model$ma)
    if (x@model$mu != 0)
        cat("\n mu:   ", x@model$mu)
    if (x@model$omega != 0)
        cat("\n omega:", x@model$omega)
    if (sum(abs(x@model$alpha)) != 0)
        cat("\n alpha:", x@model$alpha)
    if (sum(abs(x@model$gamma)) != 0)
        cat("\n gamma:", x@model$gamma)
    if (sum(abs(x@model$beta)) != 0)
        cat("\n beta: ", x@model$beta)
    if (x@model$delta != 2)
        cat("\n delta:", x@model$delta)

    # Distribution:
    cat("\nDistribution: \n ")
    cat(x@distribution)
    if (x@distribution != "norm") {
        if (x@distribution == "snorm") {
            cat("\nDistributional Parameters: \n")
            cat(" xi =", x@model$skew)
        }
        if (x@distribution == "ged" | x@distribution == "std") {
            cat("\nDistributional Parameter: \n")
            cat(" nu =", x@model$shape)
        }
        if (x@distribution == "sged" | x@distribution == "sstd") {
            cat("\nDistributional Parameters: \n")
            cat(" nu =", x@model$shape, " xi =", x@model$skew)
        }
    }

    # Seed:
    if (x@rseed != 0) {
        cat("\nRandom Seed: \n ")
        cat(x@rseed)
    }

    # Presample:
    cat("\nPresample: \n")
    n = -(length(x@presample[, 1])-1)
    time = n:0
    print(data.frame(cbind(time, x@presample)))

    # Return Value:
    invisible()
})


################################################################################

