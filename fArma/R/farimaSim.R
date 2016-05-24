
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# A copy of the GNU General Public License is available via WWW at
# http://www.gnu.org/copyleft/gpl.html.  You can also obtain it by
# writing to the Free Software Foundation, Inc., 59 Temple Place,
# Suite 330, Boston, MA  02111-1307  USA. 

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
# FUNCTIONS:            FARIMA PROCESS:
#  farimaSim             Generates FARIMA time series process
################################################################################


farimaSim = 
function(n = 1000, model = list(ar = c(0.5, -0.5), d = 0.3, ma = 0.1),
method = c("freq", "time"), ...) 
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Simulates a FARMA Time Series Process
    
    # Note:
    #   Splus-Like argument list
    
    # Example:
    #   armaSim(model = list(ar = c(0.5, -0.5), d = 0.2, ma = 0.1))
    #   armaSim(model = list(d = 0.2, ma = 0))
    #   armaSim(model = list(d = 0.2))
    
    # FUNCTION:
    
    # Settings:
    innov = NULL
    n.start = 100
    start.innov = NULL
    rand.gen = rnorm
    
    # Simulate:
    if (!is.list(model)) 
        stop("model must be list")
    if (is.null(innov)) 
        innov = rand.gen(n, ...)
    n = length(innov) 
    if (is.null(start.innov)) 
        start.innov = rand.gen(n, ...) 
    n.start = length(start.innov)

    # AR PART:
    p = length(model$ar)
    if (p == 1 && model$ar == 0) 
        p = 0
    if (p) { 
        minroots = min(Mod(polyroot(c(1, -model$ar))))
        if (minroots <= 1) stop("ar part of model is not stationary") 
    }
    
    # MA PART:
    q = length(model$ma)
    if (q == 1 && model$ma == 0) 
        q = 0
    if (n.start < p + q) 
        stop("burn-in must be as long as ar + ma")
    
    # DIFFERENCING:
    ## if (model$d < 0) stop("d must be positive ") 
    dd = length(model$d)    
    if (dd) { 
        # FRACDIFF if "dd" is a non-integer value:
        d = model$d
        if (d != round(d) ) { 
            TSMODEL = "FRACDIFF" 
        } else { 
            TSMODEL = "ARIMA" } 
    } else {
        d = 0 
        TSMODEL = "ARIMA" 
    } 
    
    # ARMA:
    if (TSMODEL == "ARIMA") {
        stop("d is a short range model")
    }
        
    if (TSMODEL == "FRACDIFF") {
        if (p == 0) model$ar = 0
        if (q == 0) model$ma = 0
        mu = 0
        # Use Fortran Routine from R's contributed fracdiff package:
        # This is a BUILTIN function ...
        x = .Fortran("fdsim", as.integer(n), as.integer(p), as.integer(q), 
            as.double(model$ar), as.double(model$ma), as.double(model$d), 
            as.double(mu), as.double(rnorm(n + q)), x = double(n + q), 
            as.double(.Machine$double.xmin), as.double(.Machine$double.xmax), 
            as.double(.Machine$double.neg.eps), as.double(.Machine$double.eps), 
            PACKAGE = "fArma")$x[1:n] 
    }
               
    # Return Value:
    ans = as.ts(x)
    attr(ans, "control") <- c(method = method, model = unlist(model))
    ans
}


################################################################################

