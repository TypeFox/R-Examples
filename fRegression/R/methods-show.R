
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


################################################################################
# S3-METHODS:           PRINT METHOD:
#  show.fREG             Prints results from a regression model fit 
################################################################################


setMethod(f = "show", signature(object = "fREG"), definition = 
    function(object) 
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Print method for Regression Modelling, an object of class "fREG"
    
    # FUNCTION:
    
    # Title:
    cat("\nTitle:\n ")
    cat(as.character(object@title), "\n")
    
    # Call:
    # cat("\nCall:\n")
    # cat(paste(deparse(object@call), sep = "\n", collapse = "\n"), 
    #    "\n", sep = "") 
        
    # Formula:
    cat("\nFormula:\n ")
    # cat(as.character(object@formula), "\n")
    print(object@formula)
    
    # Family:
    if (object@family[1] != "" && object@family[2] != "") {     
        cat("\nFamily:\n ")
        cat(as.character(object@family[1:2]), "\n") }
    
    # Digits:
    digits = max(4, getOption("digits") - 4)
        
    # Model Parameters:
    cat("\nModel Parameters:\n")        
        
        # Regression Model LM / RLM:
        if (object@method == "lm" | object@method == "rlm") {
            print.default(format(object@fit$coef, digits = digits), 
                print.gap = 2, quote = FALSE) 
        }
        
        # Regression Model GLM:
        if (object@method == "glm") {
            if (length(object@fit$coef)) {
                # if (is.character(co = object@fit$contrasts)) 
                co <- object@fit$contrasts
                if (is.character(co)) 
                    cat("  [contrasts: ", apply(cbind(names(co), co), 
                    1, paste, collapse = "="), "]")
                    # cat(":\n")
                print.default(format(object@fit$coefficients,
                    digits = digits), print.gap = 2, quote = FALSE)
            } else { 
                cat("No coefficients\n\n") 
            } 
        }
        
        # Regression Model GAM:
        if (object@method == "gam" | object@method == "am") {
            print.default(format(object@fit$coef, digits = digits), 
                print.gap = 2, quote = FALSE) 
                
        }
        
        # Regression Model PPR:
        if (object@method == "ppr") {
            cat("-- Projection Direction Vectors --\n")
            print(object@fit$alpha)
            cat("-- Coefficients of Ridge Terms --\n")
            print(object@fit$beta) 
        }                 
        
        # Regression Model POLYMARS:
        if (object@method == "polymars") {
            print(object@fit$coef) 
        }
        
        # Regression Model NNET:
        if (object@method == "nnet") {
            cat("   a ",object@fit$n[1], "-", object@fit$n[2], "-", 
                object@fit$n[3], " network", " with ", 
                length(object@fit$wts), " weights\n", sep="")
            cat("   options were -")
            tconn = diff(object@fit$nconn)
            if (tconn[length(tconn)] > object@fit$n[2]+1) 
                cat(" skip-layer connections ")
            if (object@fit$nunits > object@fit$nsunits && 
                !object@fit$softmax) 
                cat(" linear output units ")
            if (object@fit$entropy) 
                cat(" entropy fitting ")
            if (object@fit$softmax) 
                cat(" softmax modelling ")
            if (object@fit$decay[1] > 0) 
                cat(" decay=", object@fit$decay[1], sep="")
            cat("\n")
            Weights = object@fit$wts
            print(Weights) 
        } 
        
    # Residual Variance:
    # cat("\nResidual Variance:\n", var(object@fit$residuals))
    cat("\n")
        
    # Return Value:
    invisible()
})


###############################################################################


