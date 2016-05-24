
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


###############################################################################
# FUNCTION:             SIMULATION:
#  regSim                Returns a regression example data set
###############################################################################


LM3 <- 
function(n = 100, seed = 4711) 
{
    # A function implemented by Diethelm Wuertz
    
    # FUNCTION:
    
    # LM - Example Data:
    set.seed(seed)
    x1 = rnorm(n)
    x2 = rnorm(n)
    x3 = rnorm(n)
    y = 0.75 * x1 + 0.25 * x2 - 0.5 * x3
    eps = 0.1 * rnorm(n)
    y = y + eps
    data.frame(Y = y, X1 = x1, X2 = x2, X3 = x3)
}


# ------------------------------------------------------------------------------


LOGIT3 <-  
    function(n = 100, seed = 4711)
{
    # A function implemented by Diethelm Wuertz
    
    # FUNCTION:
    
    # GLM / BINOMIAL / LOGIT - Example Data:
    set.seed(seed)
    x1 = rnorm(n)
    x2 = rnorm(n)
    x3 = rnorm(n)
    eps = 0.1 * rnorm(n)
    y = 0.75 * x1 + 0.25 * x2 - 0.5 * x3 + eps
    p = 1 / ( 1 + exp(-y) )
    data.frame(Y = p, X1 = x1, X2 = x2, X3 = x3)
}


# ------------------------------------------------------------------------------


GAM3 <- 
    function(n = 100, seed = 4711)
{  
    # A function implemented by Diethelm Wuertz
    
    # FUNCTION:
    
    # GAM - Example Data:
    set.seed(seed)
    x1 = runif(n)
    x2 = runif(n)
    x3 = runif(n)
    y1 = scale(sin(2 * pi * x1))
    y2 = scale(exp(x2))
    y3 = scale(x3)
    y = scale(y1 + y2 + y3)
    eps = 0.1 * rnorm(n, sd = sd(y))
    y = y + eps
    data.frame(Y = y, X1 = x1, X2 = x2, X3 = x3)
}


# ------------------------------------------------------------------------------


regSim <- 
    function(model = "LM3", n = 100, ...)
{
    # A function implemented by Diethelm Wuertz
    
    # FUNCTION:
    
    # Simulate:
    funSim <- match.fun(model)
    ans <- funSim(n = n, ...)
    
    # Return Value:
    ans
}


###############################################################################


