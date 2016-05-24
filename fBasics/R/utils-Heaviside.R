
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
# FUNCTION:                 HEAVISIDE AND RELATED FUNCTIONS:
#  Heaviside                 Computes Heaviside unit step function
#  Sign                      Another signum function
#  Delta                     Computes delta function
#  Boxcar                    Computes boxcar function
#  Ramp                      Computes ramp function
################################################################################


Heaviside <- 
function(x, a = 0) 
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Computes the Heaviside or unit step function.
    
    # Arguments:
    #   x - a numeric vector.
    #   a - the location of the break.
    
    # Details:
    #   The Heaviside step function is 1 for x>a, 1/2 for x=a,
    #   and 0 for x<a.
    
    # Notes:
    #   Heaviside Unit Step Function and Related Functions
    #   See:  http://mathworld.wolfram.com/HeavisideStepFunction.html
    #   Note: sign(x) is part of R's base package
       
    # FUNCTION:
    
    # Compute H:
    result = (sign(x-a) + 1)/2
    
    # Return Value:
    result
}


# ------------------------------------------------------------------------------


Sign <- 
function(x, a = 0) 
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Computes the (modofied) Sign function.
    
    # Arguments:
    #   x - a numeric vector.
    #   a - the location of the break.
    
    # Details:
    #   The Sign function is 1 for x>a, 0 for x=a,
    #   and -1 for x<a.
     
    # FUNCTION:
    
    # Compute Sign:
    result = sign(x-a)
    
    # Return Value:
    result
}


# ------------------------------------------------------------------------------


Delta <- 
function(x, a = 0) 
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Computes the delta function.
    
    # Arguments:
    #   x - a numeric vector.
    #   a - the location of the break.
    
    # Details:
    #   The delta function is defined as: delta(x) = d/dx Heaviside(x-a)

    # FUNCTION:
    
    # Compute delta:
    result = 1/sign(abs(x-a))-1
    
    # Return Value:
    result
}


# ------------------------------------------------------------------------------


Boxcar <- 
function(x, a = 0.5) 
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Computes the boxcar function.
    
    # Arguments:
    #   x - a numeric vector.
    #   a - the location of the break.
    
    # Details:
    #   The boxcar function is defined as: 
    #       Pi(x) = Heaviside(x+a) - Heaviside(x-a)

    # FUNCTION:
    
    # Compute boxcar:
    result = Heaviside(x-a) - Heaviside(x+a)
    
    # Return Value:
    result
}


# ------------------------------------------------------------------------------


Ramp <- 
function(x, a = 0) 
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Computes the ramp function.
    
    # Arguments:
    #   x - a numeric vector.
    #   a - the location of the break.
    
    # Details:
    #   The ramp function is defined as: R(x)= (x-a)*Heaviside(x-a)

    # FUNCTION:
    
    # Compute ramp:
    result = (x-a) * Heaviside(x-a)
    
    # Return Value:
    result
}


################################################################################

