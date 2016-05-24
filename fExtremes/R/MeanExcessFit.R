
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
# FUNCTION:            MEAN EXCESS FUNCTION FIT:
#  normMeanExcessFit    Fits mean excesses to a normal density
#  ghMeanExcessFit      Fits mean excesses to a generalized hyperbolic density
#  hypMeanExcessFit     Fits mean excesses to a hyperbolic density
#  nigMeanExcessFit     Fits mean excesses to a normal inverse Gaussian density
################################################################################


normMeanExcessFit = 
function(x, doplot = TRUE, trace = TRUE, ...)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Fits mean excesses with a normal density   
    
    # Arguments:
    #   x -  an univariate 'timeSeries' object
    #   doplot - alogical flag. Should a mean excess plot be dispalyed?
    #   ... - optional parameters passed to the function mePlot()
    
    # FUNCTION: 
    
    # Settings:
    x = as.vector(x)
    U = mePlot(x, doplot = doplot, ...)[, 1]
    U = U[!is.na(U)]
    U = seq(min(U), max(U), length = 51)
    if(trace) print(U)
    
    # Fit Parameters:
    fit = nFit(x, doplot = FALSE, trace = FALSE)
    param = fit@fit$estimate
    
    # Compute Mean Excess Function:
    func<-function(x, u, param) {
        (x-u)*dnorm(x, param[1], param[2])}  
    Y = NULL
    for (u in U) {
        y1 = integrate(func, lower = u, upper = Inf, u = u, 
            param = param)[[1]]
        y2 = integrate(dnorm, lower = u, upper = Inf, 
            mean = param[1], 
            sd = param[2])[[1]]
        Y = c(Y, y1/y2)
    }
    
    # Plot:
    if (doplot) lines(U, Y, lwd = 2)
    
    # Result:
    result = data.frame(threshold = U, me = Y)
    attr(result, "control")<-fit 

    # Return Value:
    invisible(result)
}


# ------------------------------------------------------------------------------


ghMeanExcessFit = 
function(x, doplot = TRUE, trace = TRUE, ...)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Fits mean excesses with a hyperbolic density   
    
    # Arguments:
    #   x -  an univariate 'timeSeries' object
    #   doplot - alogical flag. Should a mean excess plot be dispalyed?
    #   ... - optional parameters passed to the function mePlot()
    
    # FUNCTION: 
    
    # Settings:
    x = as.vector(x)
    U = mePlot(x, doplot = doplot, ...)[, 1]
    U = U[!is.na(U)]
    U = seq(min(U), max(U), length = 51)
    if(trace) print(U)
    
    # Fit Parameters:
    fit = ghFit(x, doplot = FALSE, trace = FALSE)
    param = fit@fit$estimate
    
    # Compute Mean Excess Function:
    func<-function(x, u, param) {
        (x-u)*dgh(x, param[1], param[2], param[3], param[4], param[5]) }  
    Y = NULL
    for (u in U) {
        y1 = integrate(func, lower = u, upper = Inf, u = u, 
            param = param)[[1]]
        if (trace) print(c(u, y1))
        y2 = integrate(dgh, lower = u, upper = Inf, 
            alpha = param[1], 
            beta = param[2], 
            delta = param[3], 
            mu = param[4],
            lambda = param[5])[[1]]
        if (trace) print(c(u, y2))
        Y = c(Y, y1/y2)
    }
    
    # Plot:
    if (doplot) lines(U, Y, lwd = 2)
    
    # Result:
    result = data.frame(threshold = U, me = Y)
    attr(result, "control")<-fit

    # Return Value:
    invisible(result)
}


# ------------------------------------------------------------------------------


hypMeanExcessFit = 
function(x, doplot = TRUE, trace = TRUE, ...)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Fits mean excesses with a hyperbolic density   
    
    # Arguments:
    #   x -  an univariate 'timeSeries' object
    #   doplot - alogical flag. Should a mean excess plot be dispalyed?
    #   ... - optional parameters passed to the function mePlot()
    
    # FUNCTION: 
    
    # Settings:
    x = as.vector(x)
    U = mePlot(x, doplot = FALSE)[, 1]
    U = U[!is.na(U)]
    U = seq(min(U), max(U), length = 51)
    
    # Fit Parameters:
    fit = hypFit(x, doplot = FALSE, trace = FALSE)
    param = fit@fit$estimate
    
    # Compute Mean Excess Function:
    func<-function(x, u, param) {
        (x-u)*dhyp(x, param[1], param[2], param[3], param[4])}     
    Y = NULL
    for (u in U) {
        y = integrate(func, lower = u, upper = Inf, u = u, param = param)[[1]]
        Y = c(Y, y)
    }
    # Plot:
    if (doplot) lines(U, Y, lwd = 2)
    
    # Result:
    result = data.frame(threshold = U, me = Y)
    attr(result, "control")<-fit
    
    # Return Value:
    invisible(result)
}


# ------------------------------------------------------------------------------


nigMeanExcessFit = 
function(x, doplot = TRUE, trace = TRUE, ...)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Fits mean excesses with a genaralized hyperbolic density   
    
    # Arguments:
    #   x -  an univariate 'timeSeries' object
    #   doplot - alogical flag. Should a mean excess plot be dispalyed?
    #   ... - optional parameters passed to the function mePlot()
    
    # FUNCTION: 
    
    # Settings:
    x = as.vector(x)
    U = mePlot(x, doplot = doplot, ...)[, 1]
    U = U[!is.na(U)]
    U = seq(min(U), max(U), length = 51)
    if(trace) print(U)
    
    # Fit Parameters:
    fit = nigFit(x, doplot = FALSE, trace = FALSE, scale = FALSE)
    param = fit@fit$estimate
    
    # Compute Mean Excess Function:
    func<-function(x, u, param) {
        (x-u)*dnig(x, param[1], param[2], param[3], param[4]) }  
    Y = NULL
    for (u in U) {
        y1 = integrate(func, lower = u, upper = Inf, u = u, 
            param = param)[[1]]
        if (trace) print(c(u, y1))
        y2 = integrate(dnig, lower = u, upper = Inf, 
            alpha = param[1], 
            beta = param[2], 
            delta = param[3], 
            mu = param[4])[[1]]
        if (trace) print(c(u, y2))
        Y = c(Y, y1/y2)
    }
    
    # Plot:
    if (doplot) lines(U, Y, lwd = 2)
    
    # Result:
    result = data.frame(threshold = U, me = Y)
    attr(result, "control")<-fit

    # Return Value:
    invisible(result)
}


# ------------------------------------------------------------------------------


ghtMeanExcessFit = 
function(x, doplot = TRUE, trace = TRUE, ...)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Fits mean excesses with a genaralized hyperbolic density   
    
    # Arguments:
    #   x -  an univariate 'timeSeries' object
    #   doplot - alogical flag. Should a mean excess plot be dispalyed?
    #   ... - optional parameters passed to the function mePlot()
    
    # FUNCTION: 
    
    # Settings:
    x = as.vector(x)
    U = mePlot(x, doplot = doplot, ...)[, 1]
    U = U[!is.na(U)]
    U = seq(min(U), max(U), length = 51)
    if(trace) print(U)
    
    # Fit Parameters:
    fit = ghtFit(x, doplot = FALSE, trace = FALSE, scale = FALSE)
    param = fit@fit$estimate
    
    # Compute Mean Excess Function:
    func<-function(x, u, param) {
        (x-u) * dght(x, param[1], param[2], param[3], param[4]) }  
    Y = NULL
    for (u in U) {
        y1 = integrate(func, lower = u, upper = Inf, u = u, 
            param = param)[[1]]
        if (trace) print(c(u, y1))
        y2 = integrate(dght, lower = u, upper = Inf, 
            beta = param[1], 
            delta = param[2], 
            mu = param[3], 
            nu = param[4])[[1]]
        if (trace) print(c(u, y2))
        Y = c(Y, y1/y2)
    }
    
    # Plot:
    if (doplot) lines(U, Y, lwd = 2)
    
    # Result:
    result = data.frame(threshold = U, me = Y)
    attr(result, "control")<-fit

    # Return Value:
    invisible(result)
}


################################################################################

