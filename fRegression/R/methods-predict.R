
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
# FUNCTION:             DESCRIPTION REGRESSION METHODS:
#  predict.fREG          Predicts values from a fitted regression model
################################################################################


setMethod(f = "predict", signature(object = "fREG"), definition =  
    function(object, newdata, se.fit = FALSE, type = "response", ...)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Predict method for Regression Modelling, an object of class "fREG"
    
    # FUNCTION:
    
    # Fit:
    fit <- object@fit
      
    # Data as data.frame:
    if (missing(newdata)) newdata <- object@data$data
     
    # Predict:
    if (object@method == "nnet" & type == "response") type = "raw"
    ans <- .predict(object = fit, newdata = newdata, se.fit = se.fit, 
        type = type, ...) 
    
    # Make the output from 'predict' unique:
    if (se.fit) {
        if (!is.list(ans)) {
            if (is.matrix(ans)) ans = as.vector(ans)
            names(ans) = rownames(newdata) 
            ans = list(fit = ans, se.fit = NA*ans) 
        } else {
            ans = ans[1:2]
        }
    } else {
        if (is.matrix(ans)) ans = as.vector(ans)
        names(ans) = rownames(newdata) 
    }
            
    # Return Value:
    ans
})


# ------------------------------------------------------------------------------


# Note, in the following "object" concerns to the slot @fit:


.predict.lm <- function(...) stats::predict.lm(...)
    # <- function (object, newdata, se.fit = FALSE, scale = NULL, df = Inf, 
    #    interval = c("none", "confidence", "prediction"), level = 0.95, 
    #    type = c("response", "terms"), terms = NULL, na.action = na.pass, 
    #    pred.var = res.var/weights, weights = 1, ...) 

.predict.rlm <- function(...) stats::predict.lm(...)
    #

.predict.glm <- function(...) stats::predict.glm(...)
    # <- function (object, newdata = NULL, type = c("link", "response", 
    #    "terms"), se.fit = FALSE, dispersion = NULL, terms = NULL, 
    #    na.action = na.pass, ...) 

.predict.gam <- function(...) mgcv::predict.gam(...)
    # <- function (object, newdata, type = "link", se.fit = FALSE, terms = NULL, 
    #    block.size = 1000, newdata.guaranteed = FALSE, na.action = na.pass, 
    #    ...)  

.predict.ppr <- function(object, ...) { stats::predict(object, ...) }
    # <- function(object, newdata, ...)

##.predict.nnet <- function(object, ...) { nnet::predict(object, ...) }
    # <- function(object, newdata, type=c("raw","class"), ...)

##.predict.polspline <- function(object, ...) { polspline::predict(object, ...) }
    # ---- can be found in polymars.R
    # <- function(object, newdata, se.fit = FALSE, type = "response", ...)


################################################################################

