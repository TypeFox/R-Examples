
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
# You should have received A copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA


################################################################################


.distFitPlot <- 
function(fit, x, FUN = "dnig", main = "Parameter Estimation", 
    span = "auto", add = FALSE, ...)
{
    x = as.vector(x)
    if (span == "auto") span = seq(min(x), max(x), length = 501)
    z = density(x, n = 100, ...)
    x = z$x[z$y > 0]
    y = z$y[z$y > 0]
    
    # The Density function must accept multiple parameters 
    #   from the first parameter
    dFun = match.fun(FUN)
    y.points = dnig(span, fit$par)
    ylim = log(c(min(y.points), max(y.points)))
    if (add) {
        lines(x = span, y = log(y.points), col = "steelblue")
    } else {
        plot(x, log(y), xlim = c(span[1], span[length(span)]),
            ylim = ylim, type = "p", xlab = "x", ylab = "log f(x)", ...)
        title(main)
        lines(x = span, y = log(y.points), col = "steelblue")
    }
}
        
      
################################################################################
  
        