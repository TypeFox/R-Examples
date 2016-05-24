
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
# FUNCTION:                 DESCRIPION:
#  tsTest                    Time Series Test Suite
# FUNCTION:                 DEPENDENCY TEST:
#  bdsTest                   Brock-Dechert-Scheinkman test for iid series
# FUNCTION:                 NONLINEARITY TESTS:
#  wnnTest                   White Neural Network Test for Nonlinearity
#  tnnTest                   Teraesvirta Neural Network Test for Nonlinearity
################################################################################


test.tsSuite = 
function()
{  
    # NA

    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.bdsTest = 
function()
{    
    # iid example:
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")
    x = rnorm(100)
    plot(x, type = "l", col = "steelblue")
    test = bdsTest(x)
    print(test)
    p.value = as.vector(test@test$p.value)
    # Is each of the 8 p.values greater 0.1?
    checkEqualsNumeric(sum(p.value > 0.1), 8)
    
    # Not identically distributed:
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")
    x = c(rnorm(50), runif(50))
    test = bdsTest(x)
    print(test)
    p.value = as.vector(test@test$p.value)
    # Is each of the 8 p.values smaller 1e-3?
    checkEqualsNumeric(sum(p.value < 1e-3), 8)
    
    # Not independent:
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")
    n = 500
    x = rep(0, times = n)
    for(i in (2:n)) x[i] = 0.4*x[i-1] + tanh(x[i-1]) + rnorm(1, sd = 0.5)
    plot(x, type = "l", col = "steelblue")
    test = bdsTest(x)
    print(test)
    p.value = as.vector(test@test$p.value)
    # Is each of the 8 p.values smaller 1e-6?
    checkEqualsNumeric(sum(p.value < 1e-6), 8)
    

    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.wnnTest = 
function()
{    
    # White NN Test:

    # See tseries Package:
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")
    x = runif(1000, -1, 1)   
    plot(x, type = "l", col = "steelblue")
    test = wnnTest(x)
    print(test)
    p.value = as.vector(test@test$p.value)
    # Is each of the two p.values greater 0.5?
    checkTrue(as.logical(mean(p.value > 0.5)))

    ## Generate time series which is nonlinear in ``mean''
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")
    n = 1000
    x = rep(0, times = n)
    for(i in (2:n)) x[i] <- 0.4*x[i-1] + tanh(x[i-1]) + rnorm(1, sd = 0.5)
    plot(x, type = "l", col = "steelblue")
    test = wnnTest(x)
    print(test)
    p.value = as.vector(test@test$p.value)
    # Is each of the two p.values smaller than 1e-4?
    checkTrue(as.logical(mean(p.value < 1e-4)))
    
    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.tnnTest = 
function()
{   
    # Teraesvirta NN Test:

    # See example from tseries Package:
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")
    x = runif(1000, -1, 1)   
    plot(x, type = "l", col = "steelblue")
    test = tnnTest(x)
    print(test)
    p.value = as.vector(test@test$p.value)
    # Is each of the two p.values greater 0.5?
    checkTrue(as.logical(mean(p.value > 0.5)))
   
    ## Generate time series which is nonlinear in ``mean''
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")
    n = 1000
    x = rep(0, times = n)
    for(i in (2:n)) x[i] <- 0.4*x[i-1] + tanh(x[i-1]) + rnorm(1, sd = 0.5)
    plot(x, type = "l", col = "steelblue")
    test = tnnTest(x)
    print(test)
    p.value = as.vector(test@test$p.value)
    # Is each of the two p.values smaller than 1e-4?
    checkTrue(as.logical(mean(p.value < 1e-4)))
    
    # Return Value:
    return()    
}


################################################################################
    
