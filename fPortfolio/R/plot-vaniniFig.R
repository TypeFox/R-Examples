
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR Description. See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA 02111-1307 USA


################################################################################
# FUNCTION:                     DESCRIPTION:
#  .vaniniFig                    Creates Vinini's Figure in Portfolio eBook
################################################################################


.vaniniFig <- 
  function()
  {
    # A function implemented by Diethelm Wuertz
    
    # FUNCTION:
    
    # Load Dataset
    dataSet <- data("LPP2005.RET", package="fPortfolio", envir=environment())
    LPP2005.RET <- get(dataSet, envir=environment())
    
    # Example Data:
    data = 100 * LPP2005.RET[, 1:6]
    mu = colMeans(data)
    Sigma = cov(data)
    one = rep(1, times = 6)
    
    # Short Selling Solution :
    invSigma = solve(Sigma)
    A = (mu  %*% invSigma %*% one)[[1,1]]     
    C = (one %*% invSigma %*% one)[[1,1]]  
    B = (mu  %*% invSigma %*%  mu)[[1,1]]  
    E = (one %*% invSigma %*%  mu)[[1,1]]
    D = B*C - A*A
    
    # Minimum Variance Point:
    xMV = 1/sqrt(C)
    yMV = A/C
    
    # Frontier Points:
    x = seq(xMV, 4*xMV, length = 500)
    a = C
    b = -2 * A
    c = B - D * x^2
    yp = (-b + sqrt(b^2 - 4 * a * c))/(2 * a)
    ym = (-b - sqrt(b^2 - 4 * a * c))/(2 * a)
    
    # Asymptotic Slopes:
    slope = sqrt(D/C)
    intercept = A/C
    
    # Tangency Line:
    x.tg = sqrt(B)/E
    y.tg = B/E
    slope.tg = y.tg/x.tg
    
    # Plot:
    plot(x,  yp, type = "l", xlim = c(0, max(x)), ylim = c(-0.08, 0.08),
         axes = FALSE, xlab = "Covariance Risk", ylab = "Mean Retun")
    lines(x, ym)
    points(xMV, yMV, col = "orange", cex = 2, pch = 19)
    abline(intercept, slope, col = "blue", lty = 2)
    abline(intercept, -slope, col = "blue", lty = 2)
    abline(h = yMV, col = "grey", lty = 3)
    abline(v = 0, col = "grey", lty = 3)
    points(x.tg, y.tg, col = "red", cex = 2, pch = 19)
    abline(0, slope.tg, col = "brown", lty = 2)
    
    # Return Value:
    invisible()
  }


################################################################################

