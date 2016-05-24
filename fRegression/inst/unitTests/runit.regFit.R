
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
#  regFit                Wrapper Function for Regression Models
################################################################################


test.lmFit <- 
    function()
{   
    # Simulate Artificial LM:
    x <- regSim(model = "LM3", n = 50)
    
    # Fit Parameters:
    lmfit <- regFit(Y ~ X1 + X2 + X3, data = x, use = "lm") 
    
    print(lmfit)
    summary(lmfit)
    
    # plot(lmfit)
    
    fitted(lmfit)
    slot(lmfit, "fitted")
    residuals(lmfit)
    slot(lmfit, "residuals")
    
    coef(lmfit)
    formula(lmfit)
    
    predict(lmfit)
    
    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.rlmFit <- 
    function()
{   
    # Simulate Artificial LM:
    x <- regSim(model = "LM3", n = 50)
    
    # Fit Parameters:
    rlmfit <- regFit(Y ~ X1 + X2 + X3, data = x, use = "rlm") 
    
    print(rlmfit)
    summary(rlmfit)
    # plot(rlmfit)
    
    fitted(rlmfit)
    slot(rlmfit, "fitted")
    residuals(rlmfit)
    slot(rlmfit, "residuals")
    
    coef(rlmfit)
    formula(rlmfit)
    
    predict(rlmfit)
    
    head(rlmfit@fit$model)
      
    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.glmFit <- 
    function()
{   
    # Simulate Artificial LM:
    x <- regSim(model = "LM3", n = 50)
    
    # Fit Parameters:
    glmfit <- regFit(Y ~ X1 + X2 + X3, data = x, use = "glm") 
    
    print(glmfit)
    summary(glmfit)
    # plot(glmfit)
    
    print(glmfit@fit)
    summary(glmfit@fit)
    
    fitted(glmfit)
    slot(glmfit, "fitted")
    residuals(glmfit)
    slot(glmfit, "residuals")
    
    coef(glmfit)
    formula(glmfit)
    
    predict(glmfit)
     
    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.gamFit <- 
    function()
{   
    # Simulate Artificial LM:
    x <- regSim(model = "GAM3", n = 50)
    
    # Fit Parameters:
    gamfit <- regFit(Y ~ s(X1) + s(X2) + X3, data = x, use = "gam") 
    
    print(gamfit)
    summary(gamfit)
    # plot(gamfit)
    
    print(gamfit@fit)
    summary(gamfit@fit)
    
    fitted(gamfit)
    slot(gamfit, "fitted")
    residuals(gamfit)
    slot(gamfit, "residuals")
    
    coef(gamfit)
    formula(gamfit)
    
    predict(gamfit)
    
    gamfit@fit$terms
      
    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.pprFit <- 
    function()
{   
    # Simulate Artificial LM:
    x <- regSim(model = "LM3", n = 50)
    
    # Fit Parameters:
    pprfit <- regFit(Y ~ X1 + X2 + X3, data = x, use = "ppr") 
    ppr <- ppr(Y ~ X1 + X2 + X3, data = x, nterms = 2)
    
    print(pprfit)
    summary(pprfit)
    # plot(pprfit)
    
    print(pprfit@fit)
    summary(pprfit@fit)
    
    fitted(pprfit)
    slot(pprfit, "fitted")
    residuals(pprfit)
    slot(pprfit, "residuals")
    
    coef(pprfit)
    formula(pprfit)
    
    predict(pprfit)
    
    pprfit@fit$terms
       
    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


if (FALSE) { 
test.nnetFit <- 
    function()
{   
    # Simulate Artificial LM:
    x <- regSim(model = "LM3", n = 50)
    
    # Fit Parameters:
    nnetfit <- regFit(Y ~ X1 + X2 + X3, data = x, use = "nnet") 
    
    print(nnetfit)
    summary(nnetfit)
    # plot(nnetfit)
    
    print(nnetfit@fit)
    summary(nnetfit@fit)
    
    fitted(nnetfit)
    slot(nnetfit, "fitted")
    residuals(nnetfit)
    slot(nnetfit, "residuals")
    
    coef(nnetfit)
    formula(nnetfit)
    
    predict(nnetfit)
    
    nnetfit@fit$terms
       
    # Return Value:
    return()
}
}


# ------------------------------------------------------------------------------


test.polymarsFit <- 
    function()
{   
    # Simulate Artificial LM:
    x <- regSim(model = "LM3", n = 50)
    
    # Fit Parameters:
    polymarsfit <- regFit(Y ~ X1 + X2 + X3, data = x, use = "polymars") 
    
    print(polymarsfit)
    summary(polymarsfit)
    
    fitted(polymarsfit)
    slot(polymarsfit, "fitted")
    residuals(polymarsfit)
    slot(polymarsfit, "residuals")
    
    coef(polymarsfit)
    formula(polymarsfit)
    
    predict(polymarsfit)
    
    polymarsfit@fit$terms
    
    # Return Value:
    return()
}


################################################################################


