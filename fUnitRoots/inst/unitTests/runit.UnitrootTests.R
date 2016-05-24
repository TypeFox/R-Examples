
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
# FUNCTION:                ADF TESTS:
#  adfTest                  ADF unit root test using Banarjee's test statistics
#  unitrootTest             ADF unit root test using McKinnon's test statistics
################################################################################


pvalue = 
function(object) 
{   # A function implemented by Diethelm Wuertz

    # FUNCTION:
    
    # Return Value: 
    UseMethod("pvalue") 
}


# ------------------------------------------------------------------------------


pvalue.fHTEST = 
function(object) 
{   # A function implemented by Diethelm Wuertz

    # FUNCTION:
    
    # p Value:
    pValue = object@test$p.value 
    pNames = names(pValue)
    for (i in 1:length(pValue)) {
        if (is.null(pNames[i]) || pNames[i] == "") pNames[i] = "p.value"
    }
    names(pValue)<-pNames
    
    # Return Value:
    pValue        
}


################################################################################


test.adfTest = 
function()
{  
    # A time series which contains no unit-root:
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")
    x = rnorm(1000)  
    
    # A time series which contains a unit-root:
    y = cumsum(c(0, x))
    
    # Test x:
    A = adfTest(x)
    print(A)
    target = round(pvalue(A), 3)
    print(target)
    current = 0.010
    checkEqualsNumeric(target, current)
    
    # Test x:
    B = adfTest(y)
    print(B)
    target = round(pvalue(B), 3)
    print(target)
    current = 0.018
    checkEqualsNumeric(target, current)
   
    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------
 
    
test.unitrootTest = 
function()
{     
    # A time series which contains no unit-root:
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")
    x = rnorm(1000)  
    
    # A time series which contains a unit-root:
    y = cumsum(c(0, x))
    
    # Test x:
    A = unitrootTest(x)
    print(A)
    target = pvalue(A)[1]
    print(target)
    current = 0.001
    checkTrue(target < current)
    target = pvalue(A)[2]
    print(target)
    current = 0.001
    checkTrue(target < current)
    
    # Test x:
    B = unitrootTest(y)    
    print(B)
    target = pvalue(B)[1]
    print(target)
    current = 0.01
    checkTrue(target > current)
    target = pvalue(B)[2]
    print(target)
    current = 0.01
    checkTrue(target > current)
   
    # Return Value:
    return()    
}


################################################################################
    
