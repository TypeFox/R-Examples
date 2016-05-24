
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


###############################################################################


test.lmCoef <- 
    function()
{   
    # Simulate Artificial LM:
    x = regSim(model = "LM3", n = 50)
    
    # Convert to a timeSeries Object with Dummy Dates
    x = as.timeSeries(x)

    # Fit Parameters:
    fit = regFit(Y ~ X1 + X2 + X3, data = x, use = "lm") 
    fit
    
    # Extract Fitted values:
    head(slot(fit, "fitted"))
    val = fitted(fit)
    head(val)
    class(val)
    
    # Extract Residuals:
    head(slot(fit, "residuals"))
    val = residuals(fit)
    head(val)
    class(val)
    
    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.rlmCoef <- 
    function()
{   
    # Simulate Artificial LM:
    x = regSim(model = "LM3", n = 50)

    # Convert to a timeSeries Object with Dummy Dates
    x = as.timeSeries(x)
    
    # Fit Parameters:
    fit = regFit(Y ~ X1 + X2 + X3, data = x, use = "rlm") 
    fit
    
    # Extract Fitted values:
    head(slot(fit, "fitted"))
    val = fitted(fit)
    head(val)
    class(val)
    
    # Extract Residuals:
    head(slot(fit, "residuals"))
    val = residuals(fit)
    head(val)
    class(val)
    
    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.amCoef <- 
    function()
{   
    # Simulate Artificial LM:
    x = regSim(model = "GAM3", n = 50)

    # Convert to a timeSeries Object with Dummy Dates
    x = as.timeSeries(x)

    # Fit Parameters:
    fit = regFit(Y ~ X1 + X2 + X3, data = x, use = "gam") 
    fit
    
    # Extract Fitted values:
    head(slot(fit, "fitted"))
    val = fitted(fit)
    head(val)
    class(val)
    
    # Extract Residuals:
    head(slot(fit, "residuals"))
    val = residuals(fit)
    head(val)
    class(val)
    
    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.pprCoef <- 
    function()
{   
    # Simulate Artificial LM:
    x = regSim(model = "LM3", n = 50)

    # Convert to a timeSeries Object with Dummy Dates
    x = as.timeSeries(x)

    # Fit Parameters:
    fit = regFit(Y ~ X1 + X2 + X3, data = x, use = "ppr") 
    fit
    
    # Extract Fitted values:
    head(slot(fit, "fitted"))
    val = fitted(fit)
    head(val)
    class(val)
    
    # Extract Residuals:
    head(slot(fit, "residuals"))
    val = residuals(fit)
    head(val)
    class(val)
    
    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.nnetCoef <- 
    function()
{   
    # Simulate Artificial LM:
    x = regSim(model = "LM3", n = 50)

    # Convert to a timeSeries Object with Dummy Dates
    x = as.timeSeries(x)

    # Fit Parameters:
    fit = regFit(Y ~ X1 + X2 + X3, data = x, use = "nnet") 
    fit
    
    # Extract Fitted values:
    head(slot(fit, "fitted"))
    val = fitted(fit)
    head(val)
    class(val)
    
    # Extract Residuals:
    head(slot(fit, "residuals"))
    val = residuals(fit)
    head(val)
    class(val)
    
    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.polymarsCoef <- 
    function()
{   
    # Simulate Artificial LM:
    x = regSim(model = "LM3", n = 50)

    # Convert to a timeSeries Object with Dummy Dates
    x = as.timeSeries(x)

    # Fit Parameters:
    fit = regFit(Y ~ X1 + X2 + X3, data = x, use = "polymars") 
    fit
    
    # Extract Fitted values:
    head(slot(fit, "fitted"))
    val = fitted(fit)
    head(val)
    class(val)
    
    # Extract Residuals:
    head(slot(fit, "residuals"))
    val = residuals(fit)
    head(val)
    class(val)
    
    # Return Value:
    return()
}


################################################################################

