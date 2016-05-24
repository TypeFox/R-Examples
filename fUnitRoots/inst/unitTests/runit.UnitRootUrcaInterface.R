
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
# FUNCTION:                URCA [PFAGFF] UNIT ROOT TESTS:
#  urdfTest                 Augmented Dickey-Fuller test for unit roots
#  urersTest                Elliott-Rothenberg-Stock test for unit roots
#  urkpssTest               KPSS unit root test for stationarity
#  urppTest                 Phillips-Perron test for unit roots
#  urspTest                 Schmidt-Phillips test for unit roots
#  urzaTest                 Zivot-Andrews test for unit roots
################################################################################


require(urca)


pvalue = 
function(object) 
{   # A function implemented by Diethelm Wuertz

    # FUNCTION:
    
    # Return Value: 
    UseMethod("pvalue") 
}


statistic = 
function(object) 
{   # A function implemented by Diethelm Wuertz

    # FUNCTION:
    
    # Return Value: 
    UseMethod("statistic") 
}


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


statistic.fHTEST = 
function(object) 
{   # A function implemented by Diethelm Wuertz

    # FUNCTION:
    
    # Statistic:
    statistic = object@test$statistic 
    
    # Return Value:
    statistic
}


################################################################################


test.urdfTest = 
function()
{
    # A time series which contains no unit-root:
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")
    x = rnorm(100)  
    
    # A time series which contains a unit-root:
    y = cumsum(c(0, x))
    
    # DF Test:
    urdfTest(x, lags = 1, type = c("nc", "c", "ct"), 
        doplot = TRUE)                                                   

    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.urersTest = 
function()
{
    # A time series which contains no unit-root:
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")
    x = rnorm(100)  
    
    # A time series which contains a unit-root:
    y = cumsum(c(0, x))
            
    # ERS Test:
    urersTest(x, type = "DF-GLS", model = c("constant", "trend"),
        lag.max = 4, doplot = TRUE)
    urersTest(x, type = "P-test", model = c("constant", "trend"),
        lag.max = 4) # No Plot
        
    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.urkpssTest = 
function()
{
    # A time series which contains no unit-root:
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")
    x = rnorm(100)  
    
    # A time series which contains a unit-root:
    y = cumsum(c(0, x))
         
    # KPSS Test:   
    urkpssTest(x, type = "mu", lags = c("short", "long", "nil"),
        use.lag = NULL, doplot = TRUE)
    urkpssTest(x, type = "tau", lags = c("short", "long", "nil"),
        use.lag = NULL) # No Plot
      
    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.urppTest = 
function()
{
    # A time series which contains no unit-root:
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")
    x = rnorm(100)  
    
    # A time series which contains a unit-root:
    y = cumsum(c(0, x))
    
    # PP Test:  
    urppTest(x, type = "Z-alpha", model = c("constant", "trend"),
        lags = c("short", "long"), use.lag = NULL, doplot = TRUE)
    urppTest(x, type = "Z-tau", model = c("constant", "trend"),
        lags = c("short", "long"), use.lag = NULL, doplot = TRUE)
        
    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.urspTest = 
function()
{
    # A time series which contains no unit-root:
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")
    x = rnorm(100)  
    
    # A time series which contains a unit-root:
    y = cumsum(c(0, x))
    
    # SP Test:  
    urspTest(x, type = "tau", pol.deg = c(1, 2, 3, 4),
        signif = c(0.01, 0.05, 0.1), doplot = TRUE)
    urspTest(x, type = "rho", pol.deg = c(1, 2, 3, 4),
        signif = c(0.01, 0.05, 0.1), doplot = TRUE)
        
    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.urzaTest = 
function()
{        
    # A time series which contains no unit-root:
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")
    x = rnorm(100)  
    
    # A time series which contains a unit-root:
    y = cumsum(c(0, x))
    
    # ZA Test: 
    urzaTest(x, model = "intercept", lag = 2, 
        doplot = TRUE) 
    urzaTest(x, model = "trend", lag = 2, 
        doplot = TRUE) 
    urzaTest(x, model = "both", lag = 2, 
        doplot = TRUE)   

    # Return Value:
    return()    
}


################################################################################
    
