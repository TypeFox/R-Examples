
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


################################################################################
# S3-METHODS:           PLOT METHOD:    
#  plot.fREG             Plots fit and diagnostics for a regression model
#  .plot.lm               Linear Regression Model internal plot        
#  .plot.rlm              Robust Linear Regression Model internal plot
#  .plot.glm              Generalized Linear Model internal plot
#  .plot.gam              Generalized Additive Model internal plot
#  .plot.nnet             Feedforward Neural Network Model internal plot
#  .plot.ppr              Projection Pursuit Regression Model internal plot
#  .plot.polymars         Polytochomous MARS Model internal plot
# PLOTS:                DESCRIPTION:
#  .interactiveRegPlot
#  .multRegPlot
################################################################################


setMethod(f = "plot", signature(x = "fREG", y = "missing"), definition = 
    function(x, which = "ask", ...)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Plot method for an object of class 'fGARCH'
    
    # Note:
    #   This method can also be used for plotting graphs fitted by 
    #   the function 'garch' from the contributed R package 'tseries'.
    
    # FUNCTION:
        
    # Plot:
    .plot(x@fit, which = which, ...)
            
    # Return Value:
    invisible(x)
})
    
    
# ------------------------------------------------------------------------------


.plot.common <- 
    function(x, which = "ask", ...)
{
    # A function implemented by Diethelm Wuertz
    
    # Notes:
    #   1. Responses + Fitted Values Plot:
    #   2. Residuals Plot:
    #   3. Quantile Plot:
    #   4. Fitted Values vs. Residuals Plot:
    
    # FUNCTION:
    
    # Plot:
    .interactiveRegPlot(
        x,
        choices = c(
            "Responses + Fitted Values",
            "Residuals",
            "Normal Q-Q",
            "Residuals vs Fitted",
            "ACF of Residuals",
            "PACF of Residuals",
            "Positive Mean Excess Plot",
            "Negative Mean Excess Plot"),
        plotFUN = paste(".plot.", 1:8, sep = ""),
        which = which) 
            
    # Return Value:
    invisible(x)
}


# ------------------------------------------------------------------------------


.plot.1 <- function(x, ...) .responsesPlot(residuals(x)+fitted(x),fitted(x))
.plot.2 <- function(x, ...) .residualsPlot(residuals(x))    
.plot.3 <- function(x, ...)  qqnormPlot(residuals(x))
.plot.4 <- function(x, ...) .firePlot(fitted(x), residuals(x)) 
.plot.5 <- function(x, ...) .acfPlot(residuals(x))
.plot.6 <- function(x, ...) .pacfPlot(residuals(x))
.plot.7 <- function(x, ...) .mrlPlot(residuals(x))
.plot.8 <- function(x, ...) .mrlPlot(-residuals(x))


# ------------------------------------------------------------------------------


.plot.lm <- function(...) .plot.common(...)
.plot.rlm <- function(...) .plot.common(...)
.plot.glm <- function(...) .plot.common(...)
.plot.gam <- function(...) .plot.common(...)
.plot.ppr <- function(...) .plot.common(...)
.plot.nnet <- function(...) .plot.common(...)
.plot.polymars <- function(...) .plot.common(...)


# ------------------------------------------------------------------------------


.interactiveRegPlot <- 
    function(x, choices = paste("Plot", 1:19), 
    plotFUN = paste("plot.", 1:19, sep = ""), which = "all", ...)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Interactive plot method.
    
    # Arguments:
    #   x - an object to be plotted
    #   choices - the character string for the choice menu
    #   plotFUN - the names of the plot functions
    #   which - plot selection, which graph should be 
    #     displayed. If a character string named "ask" the 
    #     user is interactively asked which to plot, if
    #     a logical vector of length N, those plots which
    #     are set "TRUE" are displayed, if a character string
    #     named "all" all plots are displayed.
    
    # Note:
    #   At maximum 19 plots are supported.

    # FUNCTION:
    
    # Some cecks:
    if (length(choices) != length(plotFUN)) 
        stop("Arguments choices and plotFUN must be of same length.")
    if (length(which) > length(choices)) 
        stop("Arguments which has incorrect length.")
    if (length(which) > length(plotFUN)) 
        stop("Arguments which has incorrect length.")
    if (length(choices) > 19)
        stop("Sorry, only 19 plots at max are supported.")
    
                              
    # Plot:
    if (is.numeric(which)) {
        Which = rep(FALSE, times = length(choices))
        Which[which] = TRUE
        which = Which
    }
    if (which[1] == "all") {
        which = rep(TRUE, times = length(choices))
    }
    if (which[1] == "ask") {
        .multRegPlot(x, choices, plotFUN = plotFUN, ...) 
    } else {
        for ( i in 1:length(which) ) {
            FUN = match.fun(plotFUN[i])
            if (which[i]) FUN(x) 
        } 
    }
            
    # Return Value:
    invisible(x)
}


# ------------------------------------------------------------------------------

                
.multRegPlot <- 
    function (x, choices, plotFUN, ...) 
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    
    # Arguments:
    
    # FUNCTION:
    
    # Match Functions, up to nine(teen) ...
    if (length(plotFUN) < 19) plotFUN = 
        c(plotFUN, rep(plotFUN[1], times = 19 - length(plotFUN)))
    plot.1  = match.fun(plotFUN[1]);  plot.2  = match.fun(plotFUN[2]) 
    plot.3  = match.fun(plotFUN[3]);  plot.4  = match.fun(plotFUN[4]) 
    plot.5  = match.fun(plotFUN[5]);  plot.6  = match.fun(plotFUN[6]) 
    plot.7  = match.fun(plotFUN[7]);  plot.8  = match.fun(plotFUN[8]) 
    plot.9  = match.fun(plotFUN[9]);  plot.10 = match.fun(plotFUN[10])
    plot.11 = match.fun(plotFUN[11]); plot.12 = match.fun(plotFUN[12]) 
    plot.13 = match.fun(plotFUN[13]); plot.14 = match.fun(plotFUN[14]) 
    plot.15 = match.fun(plotFUN[15]); plot.16 = match.fun(plotFUN[16]) 
    plot.17 = match.fun(plotFUN[17]); plot.18 = match.fun(plotFUN[18]) 
    plot.19 = match.fun(plotFUN[19])        
    pick = 1
    
    while (pick > 0) { 
        pick = menu (
        ### choices = paste("plot:", choices),
        choices = paste(" ", choices), 
        title = "\nMake a plot selection (or 0 to exit):")
        # up to 19 plot functions ...
        switch (pick, 
            plot.1(x),  plot.2(x),  plot.3(x),  plot.4(x),  plot.5(x), 
            plot.6(x),  plot.7(x),  plot.8(x),  plot.9(x),  plot.10(x),
            plot.11(x), plot.12(x), plot.13(x), plot.14(x), plot.15(x), 
            plot.16(x), plot.17(x), plot.18(x), plot.19(x)) 
    } 
}


################################################################################

