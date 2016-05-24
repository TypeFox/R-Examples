
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


# fEcofin-DickeyFullerPValues.R
################################################################################
# FUNCTION:            AUGMENTED DICKEY FULLER DATA TABLES:
# adfTable              Finite sample p values for the Dickey-Fuller test
# .adfPlot              Plots sample p values for the Dickey-Fuller test
# padf                  Returns probabilities  for the ADF Test given quantiles
# qadf                  Returns quantiles for the ADF Test given probabilities
################################################################################


################################################################################
# FUNCTION:            AUGMENTED DICKEY FULLER DATA TABLES:
# adfTable              Finite sample p values for the Dickey-Fuller test
# .adfPlot              Plots sample p values for the Dickey-Fuller test
# padf                  Returns probabilities  for the ADF Test given quantiles
# qadf                  Returns quantiles for the ADF Test given probabilities


adfTable =
function(trend = c("nc", "c", "ct"), statistic = c("t", "n"), includeInf = TRUE)
{   # A function implemented by Diethelm Wuertz

    # Description:     
    #   Tables critical values for augmented Dickey-Fuller test.
    
    # Note: 
    #   x=-3:0; y=0:3; z=outer(x,y,"*"); rownames(z)=x; colnames(z)=y; z
    
    # Examples:
    #   adfTable()

    # FUNCTION:
      
    # Match Arguments:
    type = trend = match.arg(trend)
    statistic = match.arg(statistic)
    
    # Tables:
    if (statistic == "t") {
        # Hamilton Table B.6 - OLS t-Statistic
        if (type == "nc") {
            table = cbind(
                c(-2.66, -2.26, -1.95, -1.60, +0.92, +1.33, +1.70, +2.16),
                c(-2.62, -2.25, -1.95, -1.61, +0.91, +1.31, +1.66, +2.08),
                c(-2.60, -2.24, -1.95, -1.61, +0.90, +1.29, +1.64, +2.03),
                c(-2.58, -2.23, -1.95, -1.62, +0.89, +1.29, +1.63, +2.01),
                c(-2.58, -2.23, -1.95, -1.62, +0.89, +1.28, +1.62, +2.00),
                c(-2.58, -2.23, -1.95, -1.62, +0.89, +1.28, +1.62, +2.00))
        } else if (type == "c") {
            table = cbind(
                c(-3.75, -3.33, -3.00, -2.63, -0.37, +0.00, +0.34, +0.72),
                c(-3.58, -3.22, -2.93, -2.60, -0.40, -0.03, +0.29, +0.66),
                c(-3.51, -3.17, -2.89, -2.58, -0.42, -0.05, +0.26, +0.63),
                c(-3.46, -3.14, -2.88, -2.57, -0.42, -0.06, +0.24, +0.62),
                c(-3.44, -3.13, -2.87, -2.57, -0.43, -0.07, +0.24, +0.61),
                c(-3.43, -3.12, -2.86, -2.57, -0.44, -0.07, +0.23, +0.60))
        } else if (type == "ct") {
            table = cbind(
                c(-4.38, -3.95, -3.60, -3.24, -1.14, -0.80, -0.50, -0.15),
                c(-4.15, -3.80, -3.50, -3.18, -1.19, -0.87, -0.58, -0.24),
                c(-4.04, -3.73, -3.45, -3.15, -1.22, -0.90, -0.62, -0.28),
                c(-3.99, -3.69, -3.43, -3.13, -1.23, -0.92, -0.64, -0.31),
                c(-3.98, -3.68, -3.42, -3.13, -1.24, -0.93, -0.65, -0.32),
                c(-3.96, -3.66, -3.41, -3.12, -1.25, -0.94, -0.66, -0.33))      
        } else {
            stop("Invalid type specified")
        }
    } else if (statistic == "z" || statistic == "n") {
        # Hamilton Table B.5 - Based on OLS Autoregressive Coefficient
        if (type == "nc") {
            table = cbind(
                c(-11.9,  -9.3,  -7.3,  -5.3, +1.01, +1.40, +1.79, +2.28),
                c(-12.9,  -9.9,  -7.7,  -5.5, +0.97, +1.35, +1.70, +2.16),
                c(-13.3, -10.2,  -7.9,  -5.6, +0.95, +1.31, +1.65, +2.09),
                c(-13.6, -10.3,  -8.0,  -5.7, +0.93, +1.28, +1.62, +2.04),
                c(-13.7, -10.4,  -8.0,  -5.7, +0.93, +1.28, +1.61, +2.04),
                c(-13.8, -10.5,  -8.1,  -5.7, +0.93, +1.28, +1.60, +2.03))
        } else if (type == "c") {
            table = cbind(
                c(-17.2, -14.6, -12.5, -10.2, -0.76, +0.01, +0.65, +1.40),
                c(-18.9, -15.7, -13.3, -10.7, -0.81, -0.07, +0.53, +1.22),
                c(-19.8, -16.3, -13.7, -11.0, -0.83, -0.10, +0.47, +1.14),
                c(-20.3, -16.6, -14.0, -11.2, -0.84, -0.12, +0.43, +1.09),
                c(-20.5, -16.8, -14.0, -11.2, -0.84, -0.13, +0.42, +1.06),
                c(-20.7, -16.9, -14.1, -11.3, -0.85, -0.13, +0.41, +1.04))
        } else if (type == "ct") {
            table = cbind(
                c(-22.5, -19.9, -17.9, -15.6, -3.66, -2.51, -1.53, -0.43),
                c(-25.7, -22.4, -19.8, -16.8, -3.71, -2.60, -1.66, -0.65),
                c(-27.4, -23.6, -20.7, -17.5, -3.74, -2.62, -1.73, -0.75),
                c(-28.4, -24.4, -21.3, -18.0, -3.75, -2.64, -1.78, -0.82),
                c(-28.9, -24.8, -21.5, -18.1, -3.76, -2.65, -1.78, -0.84),
                c(-29.5, -25.1, -21.8, -18.3, -3.77, -2.66, -1.79, -0.87))      
        } else {
            stop("Invalid type specified")
        }
    } else {
        stop("Invalid statistic specified")
    }
            
    # Transpose:
    Table = t(table)
    colnames(Table) = c("0.010", "0.025", "0.050", "0.100", "0.900", 
        "0.950", "0.975", "0.990")
    rownames(Table) = c(" 25", " 50", "100", "250", "500", "Inf")
    ans = list(
        x = as.numeric(rownames(Table)),
        y = as.numeric(colnames(Table)),
        z = Table)
    class(ans) = "gridData"
    
    # Exclude Inf:
    if (!includeInf) {
        nX = length(ans$x)
        ans$x = ans$x[-nX]
        ans$z = ans$z[-nX, ]
    }
    
    # Add Control:
    attr(ans, "control") <-
        c(table = "adf", trend = trend, statistic = statistic)
        
    # Return Value:
    ans
} 


# ------------------------------------------------------------------------------


.adfPlot =
function(trend = c("nc", "c", "ct"), statistic = c("t", "n"))
{   # A function implemented by Diethelm Wuertz

    # Match Arguments:
    trend = match.arg(trend)
    statistic = match.arg(statistic)
    
    # Load Table:
    Y = adfTable(trend, statistic, includeInf = FALSE)
    X = cbind(expand.grid(x = Y$x, y = Y$y), z = as.vector(Y$z))
    x = X[, 1] # N
    y = X[, 3] # q-Stat
    z = X[, 2] # p-Value
    
    # Interpolate:
    gridPoints = 51
    ans = linearInterp(x, y, z, 
        xo = seq(min(x)-5, max(x)+5, length = gridPoints), 
        yo = seq(min(y), max(y), length = gridPoints)) 

    # Plot:
    persp(ans, theta = 40, xlab = "N", ylab = "q", zlab = "p-Value",
        main = paste(trend, "|", statistic, "ADF"))
    
    # Return Value:
    invisible(NULL)
}



# ------------------------------------------------------------------------------


padf = 
function(q, N = Inf, trend = c("nc", "c", "ct"), statistic = c("t", "n"))
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Returns probabilities  for the ADF Test given quantiles
    
    # Arguments:
    #   q = quantiles (-2.66 < q < 2.16)
    #   N = sample sizes
    
    # Example:
    #   padf(q = -2:2, N = 1000) 
    #   padf(q = -2:2, N = 100)
    #   padf(q = -2:2, N = 20)
    #   padf(-3, 100)
    #   q=c(-2:2); for (N in c(20, 100, 1000, 5000, Inf)) print(padf(q, N)) 
    
    # FUNCTION:
    
    # Match Arguments:
    trend = match.arg(trend)
    statistic = match.arg(statistic)
    
    # Check N:
    stopifnot(length(N) == 1)
    
    # Grid Points:
    X = adfTable(trend = trend, statistic = statistic)
    if (N < 500) {
        finiteX = cbind(expand.grid(x = X$x, y = X$y), z = as.vector(X$z))
        X = finiteX[is.finite(finiteX[, 1]), ]
        y = X[, 1] # N
        z = X[, 2] # p-Value
        x = X[, 3] # q-Stat
        # Out Points:
        xo = q
        yo = rep(N, times = length(q))
        # Interpolate:
        p = linearInterpp(x, y, z, xo = xo, yo = yo)[, 3]
    } else {
        x = X$z["Inf", ]
        y = X$y
        p = approx(x, y, xout = q)$y  
    }
    
    # Result:
    names(p) = as.character(q)
    attr(p, "control") <- c(N = N)
 
    # Return Value:
    p
} 


# ------------------------------------------------------------------------------


qadf = 
function(p, N = Inf, trend = c("nc", "c", "ct"), statistic = c("t", "n")) 
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Returns quantiles for the ADF Test given probabilities
    
    # Arguments:
    #   p = probabilities (0.01 < p < 0.99)
    #   N = sample sizes
    
    # Example:
    #   qadf(p = (1:9)/10, N = 1000) 
    #   qadf(p = (1:9)/10, N = 100)
    #   qadf(p = (1:9)/10, N = 20)
    #   qadf(0.001, 100)
    #   qadf(0.001, 1000)
    
    # FUNCTION:
    
    # Match Arguments:
    trend = match.arg(trend)
    statistic = match.arg(statistic)
    
    # Check N:
    stopifnot(length(N) == 1)
    
    # Grid Points:
    X = adfTable(trend = trend, statistic = statistic)
    if (N < 500) {
        finiteX = cbind(expand.grid(x = X$x, y = X$y), z = as.vector(X$z))
        X = finiteX[is.finite(finiteX[, 1]), ]
        x = X[, 1]
        y = X[, 2]
        z = X[, 3]
        # Out Points:
        xo = rep(N, times = length(p))
        yo = p
        # Interpolate:
        q = linearInterpp(x, y, z, xo = xo, yo = yo)[, 3]
    } else {
        x = X$y
        y = X$z["Inf", ]
        q = approx(x, y, xout = p)$y    
    }
    
    # Result:
    names(q) = as.character(p)
    attr(q, "control") <- c(N = N)
    
    # Return Value:
    q
} 


################################################################################

