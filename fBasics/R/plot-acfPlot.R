
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
# FUNCTION:             DESCRIPTION:
#  acfPlot               Displays tailored autocorrelations function plot
#  pacfPlot              Displays tailored partial autocorrelation function plot
#  teffectPlot           Estimates and displays the Taylor effect
#  lacfPlot              Displays lagged autocorrelations
################################################################################


acfPlot <-
function(x, labels = TRUE, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Autocorrelations function plot

    # Arguments:
    #   x - an uni- or multivariate return series of class 'timeSeries'
    #       or any other object which can be transformed by the function
    #       'as.timeSeries()' into an object of class 'timeSeries'.
    #   labels - a logical flag, by default true. Should a default
    #       main title and labels addet to the plot?

    # FUNCTION:

    # Settings:
    if (!is.timeSeries(x)) x = as.timeSeries(x)
    Units = colnames(x)

    # Labels:
    if (labels) {
        main = ""
        xlab = "lag"
        ylab = "ACF"
    } else {
        main = xlab = ylab = ""
    }

    # ACF:
    for (i in 1:dim(x)[2]) {
        if (labels) main = Units[i]
        ans = acf(x = as.vector(x[, i]),
            main = main, xlab = xlab, ylab = ylab, ...)
    }

    # Return Value:
    invisible(ans)
}


# ------------------------------------------------------------------------------


pacfPlot <-
function(x, labels = TRUE, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Partial autocorrelation function plot

    # Arguments:
    #   x - an uni- or multivariate return series of class 'timeSeries'
    #       or any other object which can be transformed by the function
    #       'as.timeSeries()' into an object of class 'timeSeries'.
    #   labels - a logical flag, by default true. Should a default
    #       main title and labels addet to the plot?

    # FUNCTION:

    # Settings:
    if (!is.timeSeries(x)) x = as.timeSeries(x)
    Units = colnames(x)

    # Labels:
    if (labels) {
        main = ""
        xlab = "lag"
        ylab = "PACF"
    } else {
        main = xlab = ylab = ""
    }

    # Partial ACF:
    for (i in 1:dim(x)[2]) {
        if (labels) main = Units[i]
        ans = pacf(x = as.vector(x[, i]),
            main = main, xlab = xlab, ylab = ylab, ...)
    }

    # Return Value:
    invisible(ans)
}



# ------------------------------------------------------------------------------


teffectPlot <-
function(x, deltas = seq(from = 0.2, to = 3.0, by = 0.2), lag.max = 10,
    ymax = NA, standardize = TRUE, labels = TRUE, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Evaluate and Display Taylor Effect

    # Arguments:
    #   x - an uni- or multivariate return series of class 'timeSeries'
    #       or any other object which can be transformed by the function
    #       'as.timeSeries()' into an object of class 'timeSeries'.
    #   labels - a logical flag, by default true. Should a default
    #       main title and labels addet to the plot?

    # FUNCTION:

    # Settings:
    if (!is.timeSeries(x)) x = as.timeSeries(x)
    Units = colnames(x)

    # Labels:
    if (labels) {
        main = ""
        xlab = "Exponent Delta"
        ylab = "Autocorrelation"
    } else {
        main = xlab = ylab = ""
    }

    # Taylor Effect:
    ans = list()
    maxDelta = NULL
    for (i in 1:dim(x)[2]) {
        if (labels) main = Units[i]
        X = as.vector(x[, i])
        # Standardize:
        if(standardize) X = (X-mean(X))/sqrt(var(X))
            data = matrix(data = rep(0, times = lag.max*length(deltas)),
            nrow = lag.max, byrow = TRUE)
        for (id in 1:length(deltas))
            data[,id] = as.double(acf(abs(X)^deltas[id], lag.max = lag.max,
                type = "corr", plot = FALSE)$acf)[2:(lag.max+1)]
        if (is.na(ymax)) ymax = max(data)

        # Plot:
        if (labels) {
            plot(deltas, data[1,], ylim = c(0, ymax), type = "n",
                main = main, xlab = xlab, ylab = ylab, ...)
        } else {
            plot(deltas, data[1,], type = "n",
                main = main, xlab = xlab, ylab = ylab, ...)
        }
        xl = 1:length(deltas)
        for (il in 1:(lag.max)) {
            yp = max(data[il, ])
            yl = xl[data[il, ] == yp]
            lines(deltas, data[il, ], col = il)
            points(deltas[yl], yp, pch = 19)
            maxDelta = c(maxDelta, deltas[yl])
            lines (c(1, 1), c(0, ymax))
        }

        # Grid:
        if (labels) {
            mtext("Taylor Effect", side = 4, adj = 0, col = "darkgrey",
                cex = 0.7)
            grid()
        }
        ans[[i]] = data
    }
    names(ans) = Units

    # Deltas for max Peak:
    ans$maxDelta = maxDelta

    # Return Value:
    invisible(ans)
}


# ------------------------------------------------------------------------------


lacfPlot <-
function(x, n = 12, lag.max = 20, type = c("returns", "values"),
    labels = TRUE, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Computes the lagged autocorrelation function

    # Arguments:
    #   x - an uni- or multivariate return series of class 'timeSeries'
    #       or any other object which can be transformed by the function
    #       'as.timeSeries()' into an object of class 'timeSeries'.
    #   type - a character string which specifies the type of the input
    #       series, either "returns" or series "values". In the case of
    #       a return series as input, the required value series is
    #       computed by cumulating the financial returns: 'exp(colCumsums(x))'
    #   labels - a logical flag, by default true. Should a default
    #       main title and labels addet to the plot?

    # FUNCTION:

    # Settings:
    if (!is.timeSeries(x)) x = as.timeSeries(x)
    Units = colnames(x)

    # Cumulated Returns:
    type = match.arg(type)
    if (type == "values") {
        cumX = x
    } else if (type == "returns") {
        cumX = exp(colCumsums(x))
    }

    # Labels:
    if (labels) {
        main = ""
        xlab = "tau"
        ylab = "Lagged Correlation"
    } else {
        main = xlab = ylab = ""
    }

    cumRho = cumLagged = NULL
    DIM = dim(cumX)[2]
    for (i in 1:DIM) {

        # Get Data:
        x = as.vector(as.matrix(cumX)[, i])
        if (labels) main = Units[i]

        # Truncate to multiple of n:
        N = trunc(length(x)/n)
        M = length(x) - n*N
        if (M > 0) x = x[-c(1:M)]

        # One Step Volatilities:
        x.ret = c(0, diff(log(x)))
        x.mat = matrix(x.ret, byrow = TRUE, ncol = n)
        u = apply(abs(x.mat), 1, mean)

        # n-step Volatilities:
        index = n*(1:N)
        v = abs(c(0, diff(log(x[index]))))

        # Zero Tau:
        L = length(u)
        RhoZero = cor(u, v)
        # print(RhoZero)

        # Positive Tau:
        RhoPos = NULL
        for (tau in 1:lag.max) {
            X = u[-((L-tau+1):L)]
            X2 = X
            Y = v[-((L-tau+1):L)]
            Y2 = v[-(1:tau)]
            X.mean = mean(X)
            Y.mean = mean(Y)
            X1 = sum((X - X.mean)^2)
            Y1 = sum((Y - Y.mean)^2)
            XY1 = sum( (X2-X.mean)*(Y2-Y.mean) )
            rho = XY1/sqrt(X1*Y1)
            RhoPos = c(RhoPos, rho)
        }

        # Negative Tau:
        RhoNeg = NULL
        for (tau in 1:lag.max) {
            X = v[-((L-tau+1):L)]
            X2 = X
            Y = u[-((L-tau+1):L)]
            Y2 = u[-(1:tau)]
            X.mean = mean(X)
            Y.mean = mean(Y)
            X1 = sum((X - X.mean)^2)
            Y1 = sum((Y - Y.mean)^2)
            XY1 = sum( (X2-X.mean)*(Y2-Y.mean) )
            rho = XY1/sqrt(X1*Y1)
            RhoNeg = c(RhoNeg, rho)
        }

        # Correlations:
        Lagged = RhoPos - RhoNeg
        Rho = c(rev(RhoNeg), RhoZero, RhoPos)

        # Plot:
        plot(x = (-lag.max):(lag.max), y = Rho, type = "l", xlab = xlab,
            ylab = ylab, ylim = c(min(Lagged), max(Rho)),
            main = main, ...)
        Text = paste("n =", n, " |  lag =", lag.max)
        mtext(Text, side = 4, adj = 0, col ="darkgrey", cex = 0.7)
        points(-lag.max:lag.max, Rho, pch = 19, cex = 0.7)
        lines(0:lag.max, c(0, Lagged), col = "red")
        points(0:lag.max, c(0, Lagged), pch = 19, cex = 0.7, col = "red")
        abline(h = 0, col = "grey", lty = 3)
        ci = 1/sqrt(length(u))
        abline(h = +ci, col = "blue")
        abline(h = -ci, col = "blue")
        if (labels) grid()

        # Grid:
        if (labels) grid()

        cumRho = rbind(cumRho, Rho)
        cumLagged = c(cumLagged, Lagged)

    }

    # Return Value:
    invisible(list(Rho = cumRho, Lagged = cumLagged))
}


################################################################################

