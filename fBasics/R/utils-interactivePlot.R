
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
# FUNCTION:                 PLOT UTILITIES:
#  interactivePlot           Plots several graphs interactively
################################################################################


interactivePlot <-
function(x, choices = paste("Plot", 1:9),
    plotFUN = paste("plot.", 1:9, sep = ""), which = "all", ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Plot method for an object of class "template".

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
    #   At maximum 9 plots are supported.

    # FUNCTION:

    # Some cecks:
    if (length(choices) != length(plotFUN))
        stop("Arguments choices and plotFUN must be of same length.")
    if (length(which) > length(choices))
        stop("Arguments which has incorrect length.")
    if (length(which) > length(plotFUN))
        stop("Arguments which has incorrect length.")
    if (length(choices) > 9)
        stop("Sorry, only 9 plots at max are supported.")

    # Internal "askPlot" Function:
    multPlot = function(x, choices, ...)
    {
        # Selective Plot:
        selectivePlot <-
    function(x, choices, FUN, which){
            # Internal Function:
            askPlot <-
        function(x, choices, FUN) {
                # Pick and Plot:
                pick = 1
                setRmetricsOptions(n.plots = length(choices))
                while (pick > 0) { pick = menu (
                    choices = paste("plot:", choices),
                    title = "\nMake a plot selection (or 0 to exit):")
                    if (pick > 0) match.fun(FUN[pick])(x) }
            }
            if (as.character(which[1]) == "ask") {
                askPlot(x, choices = choices, FUN = FUN, ...)
            } else {
                for (i in 1:getRmetricsOptions("n.plots"))
                    if (which[i]) match.fun(FUN[i])(x)
            }
            invisible()
        }

        # match Functions, up to nine ...
        if (length(plotFUN) < 9) plotFUN =
            c(plotFUN, rep(plotFUN[1], times = 9 - length(plotFUN)))
        plot.1 = match.fun(plotFUN[1]); plot.2 = match.fun(plotFUN[2])
        plot.3 = match.fun(plotFUN[3]); plot.4 = match.fun(plotFUN[4])
        plot.5 = match.fun(plotFUN[5]); plot.6 = match.fun(plotFUN[6])
        plot.7 = match.fun(plotFUN[7]); plot.8 = match.fun(plotFUN[8])
        plot.9 = match.fun(plotFUN[9])
        pick = 1
        while (pick > 0) { pick = menu (
            ### choices = paste("plot:", choices),
            choices = paste(" ", choices),
            title = "\nMake a plot selection (or 0 to exit):")
            # up to 9 plot functions ...
            switch (pick, plot.1(x), plot.2(x), plot.3(x), plot.4(x),
                plot.5(x), plot.6(x), plot.7(x), plot.8(x), plot.9(x) )
        }
    }

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
        multPlot(x, choices, ...)
    } else {
        for ( i in 1:length(which) ) {
            FUN = match.fun(plotFUN[i])
            if (which[i]) FUN(x)
        }
    }

    # Return Value:
    invisible(x)
}


################################################################################

