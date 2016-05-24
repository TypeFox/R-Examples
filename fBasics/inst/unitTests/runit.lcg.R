
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
#   1999 - 2008, Diethelm Wuertz, Rmetrics Foundation, GPL
#   Diethelm Wuertz <wuertz@itp.phys.ethz.ch>
#   www.rmetrics.org
# for the code accessed (or partly included) from other R-ports:
#   see R's copyright and license files
# for the code accessed (or partly included) from contributed R-ports
# and other sources
#   see Rmetrics's copyright file 


################################################################################
# FUNCTION:            PORTABLE INNOVATIONS:
#  set.lcgseed          Sets initial random seed
#  get.lcgseed          Gets the current valus of the random seed
# FUNCTION:            DISTRIBUTIONS:
#  runif.lcg            Generates portable uniform linear congruational rvs
#  rnorm.lcg            Generates portable normal linear congruational rvs
#  rt.lcg               Generates portable Student-t linear congruational rvs
################################################################################


test.portableInnovations <- 
    function()
{
    #  set.lcgseed  Sets initial random seed
    #  get.lcgseed  Gets the current valus of the random seed
    
    # Set initial random seed - set.lcgseed(seed = 4711)
    set.lcgseed()
    
    # Get the current valus of the random seed - get.lcgseed()
    get.lcgseed()
    
    # Return Value:
    return()
}

    
# ------------------------------------------------------------------------------


test.randomNumbers <- 
    function()
{
    #  runif.lcg    Generates portable uniform linear congruational rvs
    #  rnorm.lcg    Generates portable normal linear congruational rvs
    #  rt.lcg       Generates portable Student-t linear congruational rvs
    
    # Uniform:
    x = runif.lcg(n = 1000, min = 0, max = 1)
    hist(x)
    
    # Normal:
    x = rnorm.lcg(n = 1000, mean = 0, sd = 1)
    hist(x)
    
    # Student-t:
    x = rt.lcg(n = 1000, df = 4)
    hist(x)
    
    # Return Value:
    return()
}


################################################################################

