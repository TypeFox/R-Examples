
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
# FUNCTION:             REGRESSION MODELLING DESCRIPTION:
#  polymars              Polymars Regression
#  .polymars             Polymars regress from package polspline
#  .polymarsDefault      Internal Function
#  .polymarsVormula      Internal Function
#  .predict.polymars     Internal Function
################################################################################


test.polymars <- 
    function()
{   
    # Simulate Artificial LM:
    x <- regSim(model = "LM3", n = 50)
    
    # Original Polymars:
    fit <- polspline::polymars(responses = x[,1], predictors = x[,2:4])
    
    # Model Fitting:
    fit$fitting
    
    # Model Produced:
    fit$model
    fit$coef
    
    # Summary:
    #   Note print.polymars = summary.polymars
    polspline::summary.polymars(fit)
    
    # Predict:
    ans <- polspline::predict.polymars(object = fit, x = x[,-1]) 
    as.vector(ans)
    as.vector(fit$fitted)
    
    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.polymarsDefault <- 
    function()
{   
    # Simulate Artificial LM:
    set.seed(4711)
    x = regSim(model = "LM3", n = 50)
    
    # Polymars Wrapper:
    fit1 = fRegression:::.polymarsDefault(responses = x[,1], predictors = x[, 2:4])
    class(fit1)
    names(fit1)
    
    # Note, this fails:
    # fit1 = .polymars(responses = x[,1], predictors = x[,2:4])
    
    # Model Fitting:
    fit1$fitting
    
    # Model Produced:
    #   fit1$model reserved for model series, use ...
    fit1$coef
    
    # Summary:
    print(fit1)
    
    # Print:
    summary(fit1)
    
    # Predict:
    ans <- polspline::predict.polymars(object = fit1, x = x[,-1]) 
    as.vector(ans)
    as.vector(fit1$fitted)
    
    # Check:
    fit1$ranges.and.medians 
       
    # Return Value:
    return()
}


# -----------------------------------------------------------------------------


test.polymarsFormula <- 
    function()
{   
    # Simulate Artificial LM:
    set.seed(4711)
    x <- regSim(model = "LM3", n = 50)
    
    # Polymars Formula Wrapper:
    fit2 <- fRegression:::.polymarsFormula(formula = Y ~ X1 + X2 + X3, data = x)
    fit2 <- fRegression:::.polymars(formula = Y ~ X1 + X2 + X3, data = x)
    class(fit2)
    names(fit2)
    
    # Model Fitting:
    fit2$fitting
    
    # Model Produced:
    #   fit$model reserved for model series, use ...
    fit2$coef
    
    # Summary:
    print(fit2)
    
    # Print:
    summary(fit2)
    
    # Predict:
    fit2$model <- fit2$coef
    ans <- polspline::predict.polymars(object = fit2, x = x[,-1]) 
    as.vector(ans)
    as.vector(fit2$fitted)
    
    # Check:
    fit2$ranges.and.medians 
    
    # Return Value:
    return()
}

# -----------------------------------------------------------------------------


test.regFit.polymars <- 
    function()
{   
    # Simulate Artificial LM:
    x <- regSim(model = "LM3", n = 50)
    
    # Polymars Formula Wrapper:
    fit <- regFit(formula = Y ~ X1 + X2 + X3, data = x, use = "polymars")
    class(fit)
    
    # Model Fitting:
    fit@fit$fitting
    
    # Model Produced:
    #   fit$model reserved for model series, use ...
    fit@fit$coef
    
    # Summary:
    print(fit)
    
    # Print:
    summary(fit)
           
    # Return Value:
    return()
}


# -----------------------------------------------------------------------------


test.regFit.polymars.methods <- 
    function()
{   
    # Simulate Artificial LM:
    x <- regSim(model = "LM3", n = 20)
    
    # Polymars Formula Wrapper:
    polymarsfit <- regFit(
      formula = Y ~ X1 + X2 + X3, data = x, use = "polymars")
   
    # Print:
    print(polymarsfit)
    
    # Summary:
    summary(polymarsfit)
    
    # Fitted Values:
    fitted(polymarsfit)
    slot(polymarsfit, "fitted")
    
    # Residuals:
    residuals(polymarsfit)
    slot(polymarsfit, "residuals")
    
    # Coefficients:
    coef(polymarsfit)
    
    # Formula
    formula(polymarsfit)
    
    # Return Value:
    return()
}


# -----------------------------------------------------------------------------


test.regFit.polymars.predict <- 
    function()
{   
    # Simulate Artificial LM:
    x <- regSim(model = "LM3", n = 50)
    
    # regFit / Polymars Formula Wrapper:
    fit <- regFit(formula = Y ~ X1 + X2 + X3, data = x, use = "polymars")
    class(fit)
    fit@fit$cmd
     
    # Predict from predict.polymars:
    object <- fit@fit
    class(object) = "polymars"
    object
    object$model = object$coef
    ans <- polspline::predict.polymars(object = object, x = x[,-1]) 
    as.vector(ans)
    as.vector(fit@fitted)
    
    # Predict from predict.fREG:
    ans <- predict(object = fit, newdata = x) 
    as.vector(ans)
    as.vector(fit@fitted)
    
    # Return Value:
    return()
}


###############################################################################


