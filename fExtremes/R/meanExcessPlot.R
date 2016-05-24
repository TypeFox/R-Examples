
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
# FUNCTION:            MEAN EXCESS FUNCTION PLOT:
#  meanExcessPlot       Plot mean excesses to a normal/nig/ght density
################################################################################


.meanExcessPlot <-
    function(x, labels = TRUE, title = FALSE, grid = TRUE,
    col = "steelblue", ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Plots and fits mean excess function

    # Arguments:

    # FUNCTION:

    # Common Range:
    DIM = NCOL(x)
    xRange = yRange = NULL
    for (i in 1:DIM) {
        xRange = c(xRange,
            range(mePlot(-scale(x[, i]), doplot = FALSE)[, 1], na.rm = TRUE))
        yRange = c(yRange,
            range(mePlot(-scale(x[, i]), doplot = FALSE)[, 2], na.rm = TRUE))
    }
    xLim = range(xRange)
    yLim = c(0, max(yRange))
    xPos = min(xLim) + 0.075*diff(xLim)
    yPos = 0.05*diff(yLim)

    # Colors:
    if (length(col) == 1) col = rep(col, times = DIM)

    # Labels:
    if (title) {
        xlab = "Threshold"
        ylab = "Mean Excess"
        main = colnames(X)
    } else {
        xlab = ylab = main = ""
    }

    # Mean Excess:
    for (i in 1:DIM)
    {
        # Scale Tail of Series:
        X = -scale(x[, i])
        if (labels) main = colnames(X)

        # Normal Fit:
        me = normMeanExcessFit(X, doplot = TRUE, trace = FALSE, lwd = 2,
            labels = FALSE, col = col[i], xlim = xLim, ylim = yLim,
            main = main, xlab = xlab, ylab = ylab, ...)

        normLLH = attr(me, "control")@fit$minimum

        if (grid) {
            grid(col = "darkgrey")
        }

        if (title) {
            mtext("Scaled Mean Excess", line = 0.5, cex = 0.7)
        }

        # Add 95% and 99% Sample Quantiles:
        abline(v = quantile(X, 0.95, type = 1), col = "darkgrey")
        abline(v = quantile(X, 0.99, type = 1), col = "darkgrey")

        # If Normality rejected, add  NIG and GH Student-t:
        test = jbTest(X)@test$p.value[3]
        nigLLH = ghtLLH = -9.99e99
        if (test == 0)
        {
            # NIG Fit:
            me = nigMeanExcessFit(X, doplot = FALSE, trace = FALSE)
            lines(me, col = "green", lwd = 2)
            nigLLH = attr(me, "control")@fit$minimum

            param = attr(me, "control")@fit$estimate
            abline(v = qnig(0.95, param[1], param[2], param[3], param[4]),
                col = "green")
            abline(v = qnig(0.99, param[1], param[2], param[3], param[4]),
                col = "green")

            # GH Student-t Fit:
            me = ghtMeanExcessFit(X, doplot = FALSE, trace = FALSE)
            lines(me, col = "red", lwd = 2)
            ghtLLH = attr(me, "control")@fit$minimum
        }

        # Finish:
        if (title) {
            LLH = c("NORM", "NIG", "GHT")
            colorsLLH = c("black", "green", "red")
            if (test == 0) {
                mText = paste(
                    "logLLH: NORM = ", signif(normLLH, 5),
                    " | NIG = ", signif(nigLLH, 5),
                    " | GHT = ", signif(ghtLLH, 5), sep = "")
                mtext(mText, side = 4, adj = 0, col = "darkgrey", cex = 0.7)
            } else {
                mText = paste(
                    "logLLH: NORM = ", signif(normLLH, 5), sep = "")
                mtext(mText, side = 4, adj = 0, col = "darkgrey", cex = 0.7)
            }
            indexLLH = which.max(c(normLLH, nigLLH, ghtLLH))
            maxLLH = LLH[indexLLH]
            colLLH = colorsLLH[indexLLH]
            text(xPos, yPos, maxLLH, col = colLLH)
        }
    }

    # Return Value:
    invisible()
}


################################################################################

