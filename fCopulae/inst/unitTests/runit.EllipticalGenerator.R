
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
# FUNCTION:                  UTILITY FUNCTIONS:
#  ellipticalList             Returns list of implemented Elliptical copulae
#  ellipticalParam            Sets default parameters for an elliptical copula
#  ellipticalRange            Returns the range of valid rho values
#  ellipticalCheck            Checks if rho is in the valid range
# FUNCTION:                  ELLIPTICAL GENERATOR AND RELATED FUNCTIONS:
#  gfunc                      Generator function for elliptical distributions
#  gfuncSlider                Slider for generator, density and probability
#  .pelliptical               Univariate elliptical distribution probability
#  .delliptical               Univariate elliptical distribution density
#  .qelliptical               Univariate elliptical distribution quantiles
################################################################################


test.ellipticalList = 
function()
{
    # Arguments ?
    args(ellipticalList)
    
    # List:
    target = c("norm", "cauchy", "t", "logistic", "laplace", "kotz", "epower")
    current = ellipticalList()
    print(current)
    
    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.ellipticalRange = 
function()
{
    # Arguments ?
    args(ellipticalRange)
    
    # Range:
    for (type in ellipticalList()) {
        cat("\n")
        print(ellipticalRange(type))
    }
    
    # Return Value:
    return()    
}



# ------------------------------------------------------------------------------


test.ellipticalParam = 
function()
{
    # Arguments ?
    args(ellipticalParam)
    
    # Parameters:
    for (type in ellipticalList()) {
        cat("\n")
        print(unlist(ellipticalParam(type)))
    }
    
    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.ellipticalCheck = 
function()
{
    # Arguments ?
    args(ellipticalCheck)
    # ellipticalCheck(rho = 0.75, param = NULL, type = ellipticalList()) 

    
    # Range:
    for (type in ellipticalList()) {
        cat("\n")
        param = ellipticalParam(type)$param
        rho = param[1]
        # Returns NULL if OK
        print(ellipticalCheck(rho, param[-1], type))        
    }

    # Return Value:
    return()    
}


################################################################################


test.gfunc = 
function()
{
    # Arguments ?
    args(gfunc)
    # gfunc(x, param = NULL, type = ellipticalList()) 


    # Call Generator Function - Missing x:
    gfunc(type = "norm")
    gfunc(type = "cauchy")
    gfunc(type = "t")
    gfunc(type = "t", param = 2)
    gfunc(type = "logistic")
    gfunc(type = "laplace")
    gfunc(type = "kotz")
    gfunc(type = "kotz", param = 2)
    gfunc(type = "epower")  
    gfunc(type = "epower", param = c(2, 1))
      
    # Call Generator Function - With specified x:
    gfunc(x = 0:10, type = "norm")
    gfunc(x = 0:10, type = "cauchy")
    gfunc(x = 0:10, type = "t") 
    gfunc(x = 0:10, type = "logistic")
    gfunc(x = 0:10, type = "laplace")
    gfunc(x = 0:10, type = "kotz") 
    gfunc(x = 0:10, type = "epower") 
    
    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.gfuncSlider = 
function()
{
    # Try Slider:
    # gfuncSlider()
    NA
    
    # Return Value:
    return()    
}

    
# ------------------------------------------------------------------------------


test.pelliptical = 
function()
{  
    # Probability:
    q = (-1000:1000)/2000
    S = NULL
    
    s = Sys.time()
    fCopulae:::.pelliptical(q = q, param = NULL, type = "norm")
    S = c(S, as.integer(Sys.time() - s))
     
    s = Sys.time()
    fCopulae:::.pelliptical(q = q, param = NULL, type = "cauchy") 
    S = c(S, as.integer(Sys.time() - s))
    
    s = Sys.time()
    fCopulae:::.pelliptical(q = q, param = 2, type = "t") 
    S = c(S, as.integer(Sys.time() - s))
    
    s = Sys.time()
    fCopulae:::.pelliptical(q = q, param = NULL, type = "logistic") 
    S = c(S, as.integer(Sys.time() - s))
    
    s = Sys.time()
    fCopulae:::.pelliptical(q = q, param = NULL, type = "laplace")
    S = c(S, as.integer(Sys.time() - s))
    
    s = Sys.time()
    fCopulae:::.pelliptical(q = q, param = c(r = 1), type = "kotz")  
    S = c(S, as.integer(Sys.time() - s))
    
    s = Sys.time()
    fCopulae:::.pelliptical(q = q, param = c(r = 1, s = 1), type = "epower")
    S = c(S, as.integer(Sys.time() - s))
    
    # Return Value:
    return() 
}


# ------------------------------------------------------------------------------


test.delliptical = 
function()
{   
    # Probability:
    N = 100
    x = (-1999:1999)/N
    d = fCopulae:::.delliptical(x = x, param = NULL, type = "norm")
        sum(d)/N
    d = fCopulae:::.delliptical(x = x, param = NULL, type = "cauchy")
        sum(d)/N
    d = fCopulae:::.delliptical(x = x, param = NULL, type = "t")
        sum(d)/N
    d = fCopulae:::.delliptical(x = x, param = NULL, type = "logistic")
        sum(d)/N
    d = fCopulae:::.delliptical(x = x, param = NULL, type = "laplace")
        sum(d)/N
    d = fCopulae:::.delliptical(x = x, param = NULL, type = "kotz")
        sum(d)/N
    d = fCopulae:::.delliptical(x = x, param = NULL, type = "epower")
        sum(d)/N
        
    # Non-default Parameters:
    d = fCopulae:::.delliptical(x = (-100:100)/10, param = 1, type = "kotz")
        sum(d)/N
    d = fCopulae:::.delliptical(x = (-100:100)/10, param = 1/2, type = "kotz")
        sum(d)/N
        
    # Return Value:
    return() 
}


# ------------------------------------------------------------------------------


test.qelliptical = 
function()
{   
    # Probability:
    p = (0:10)/10
    fCopulae:::.qelliptical(p = p, param = NULL, type = "norm") 
    fCopulae:::.qelliptical(p = p, param = NULL, type = "cauchy") 
    fCopulae:::.qelliptical(p = p, param = 2, type = "t") 
    fCopulae:::.qelliptical(p = p, param = NULL, type = "logistic") 
    fCopulae:::.qelliptical(p = p, param = NULL, type = "laplace")
    fCopulae:::.qelliptical(p = p, param = c(r = 1), type = "kotz")  
    fCopulae:::.qelliptical(p = p, param = c(r = 1, s = 1), type = "epower") 
    
    # Return Value:
    return() 
}


################################################################################

