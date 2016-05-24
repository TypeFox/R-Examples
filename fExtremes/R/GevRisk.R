
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
# FUNCTION:             ADDITIONAL FUNCTIONS:
#  gevrlevelPlot         Calculates Return Levels Based on GEV Fit
#  .gevrlevelLLH         Computes log-likelihood function for gevrlevelPlot
################################################################################


gevrlevelPlot =
function(object, kBlocks = 20,  ci = c(0.90, 0.95, 0.99), 
plottype = c("plot", "add"), labels = TRUE,...)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Calculates Return Levels Based on GEV Fit
    
    # Arguments:
    #   object - an object of class "fGEVFIT" as returned by the 
    #       function gevFit().
    #   kBlocks - specifies the particular return level to be 
    #       estimated; default set arbitrarily to 20

    # Note:
    #   Partial copy from R package evir
    
    # Examples:
    #   ans = gevFit(gevSim(), type = "mle", gumbel = FALSE)
    #   ans = gevrlevelPlot(ans); ans@fit$rlevel
    #   ans = gevFit(.gumbelSim(), type = "mle", gumbel = TRUE)
    #   ans = gevrlevelPlot(ans); ans@fit$rlevel
    #
    #   BMW annual (12 month) Return Level: 
    #   ans = gevFit(as.timeSeries(data(bmwRet)), "m"); gevrlevelPlot(ans, 12)
    
    # FUNCTION:
    
    # Check:
    stopifnot(object@method[1] == "gev")
    stopifnot(object@method[2] == "mle")
    stopifnot(kBlocks > 1)
    stopifnot(max(ci) < 1)
    stopifnot(min(ci) > 0)
    
    # Settings:
    out = object@fit
    conf = ci[1]
    plottype = plottype[1]
    
    # Data:
    par.ests = out$par.ests
    mu = par.ests["mu"]
    beta = par.ests["beta"]
    xi = par.ests["xi"]
    pp = 1/kBlocks
    v = qgev((1 - pp), xi, mu, beta)
    if (plottype[1] == "add") abline(h = v)
    data = out$data
    overallmax = out$llh # DW: out$nllh.final
    
    beta0 = sqrt(6 * var(data))/pi
    xi0 = 0.01
    theta = c(xi0, beta0)
    
    # Return Levels:
    parmax = NULL
    rl = v * c(0.5, 0.6, 0.7, 0.8, 0.85, 0.9, 0.95, 1, 1.1, 1.2,
        1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3, 3.25, 3.5, 4.5)
    for (i in 1:length(rl)) {
        fit = optim(theta, .gevrlevelLLH, tmp = data, pp = pp, rli = rl[i])
        parmax = rbind(parmax, fit$value)
    }
    parmax = -parmax
    overallmax = -overallmax
    crit = overallmax - qchisq(0.9999, 1)/2
    cond = parmax > crit
    rl = rl[cond]
    parmax = parmax[cond]
    smth = spline(rl, parmax, n = 200)
    aalpha = qchisq(conf[1], 1)
    
    # Labels:
    if (labels) {
        main = paste(kBlocks, "Blocks Return Level")
        xlab = "rl"
        ylab = "parmax"
    } else {
        main = xlab = ylab = ""
    }
    
    # Plot ?
    if (plottype[1] == "plot") {
        plot(rl, parmax, type = "p", pch = 19, col = "steelblue",
            main = main, xlab = xlab, ylab = ylab, ...)
        h = overallmax - aalpha/2
        abline(h = h, lty = 3, col = "brown")
        abline(v = v, lty = 3, col = "brown")
        lines(smth, ...)
        if (labels) {
            ciText = paste(as.character(100*conf[1]), "%", sep = "")
            span = 0.025*abs(max(parmax)-min(parmax))
            text(max(rl), h+span, ciText)
        }
        if (length(ci) > 1) {
            for ( i in 2:length(ci) ) {
                gevrlevelPlot(object = object, kBlocks = kBlocks, 
                    ci = ci[i], plottype = c("nextconf"), 
                    labels = labels, ...)
            }
        }
    }
    
    # Internally used to add furter confidence level lines ...
    if (plottype[1] == "nextconf") {
        h = overallmax - aalpha/2
        abline(h = h, lty = 3, col = "brown")
        abline(v = v, lty = 3, col = "brown")
        lines(smth, ...)
        if (labels) {
            ciText = paste(as.character(100*conf[1]), "%", sep = "")
            span = 0.025*abs(max(parmax)-min(parmax))
            text(max(rl), h+span, ciText)
        }
    }
    
    # Or Add ?
    ind = smth$y > overallmax - aalpha/2
    ci = range(smth$x[ind])
    if (plottype[1] == "add") {
        abline(v = ci[1], lty = 2, col = "orange")
        abline(v = ci[2], lty = 2, col = "orange")
    }
    
    # Result:
    ans = as.numeric(c(ci[1], v, ci[2]))
    ans = data.frame(cbind(min = ans[1], v = ans[2], max = ans[3], 
        kBlocks = kBlocks), row.names = "GEV Return Level")
    
    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


.gevrlevelLLH = 
function(theta, tmp, pp, rli)
{   # A copy from evir

    # Description:
    #   Computes log-likelihood function for gevrlevelPlot
    
    # Arguments:
    
    # FUNCTION:
    
    # LLH:
    mu = rli + (theta[2]*(1-(-log(1-pp))^(-theta[1])))/theta[1]
    y = 1 + (theta[1]*(tmp-mu))/theta[2]
    if ((theta[2] < 0) | (min(y) < 0)) {
        ans = NA
    } else {
        term1 = length(tmp) * log(theta[2])
        term2 = sum((1 + 1/theta[1]) * log(y))
        term3 = sum(y^(-1/theta[1]))
        ans = term1 + term2 + term3
    }
    
    # Return Value:
    ans
}


################################################################################

