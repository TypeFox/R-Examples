
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
# FUNCTION:                 DESCRIPTION:
# .jbTable                   Finite sample p values for the Jarque Bera test
# .jbPlot                    Plots probabilities
# .pjb                       Returns probabilities for JB given quantiles
# .qjb                       Returns quantiles for JB given probabilities
################################################################################


.jbTable <- 
function(type = c("LM", "ALM"), size = c("mini", "small", "all"))
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Finite sample p values for the Jarque Bera test

    # Details:
    #   The function jbTable() returns a data.frame with rows denoting
    #   size and columns denoting probabilities 0 < p < 1.

    # Note:
    #   x=-3:0; y=0:3; z=outer(x,y,"*"); rownames(z)=x; colnames(z)=y; z

    # Example:
    #   .jbTable()

    # FUNCTION:

    # Check Arguments:
    type = match.arg(type)
    size = match.arg(size)

    # Create Table:
    if (type == "LM") {
        Table = t(.jbLM())
    } else if (type == "ALM") {
        Table = t(.jbALM())
    }

    # Select Size:
    if (size == "small") {
        Table = t(Table)
        n = dim(Table)[1]
        Table = Table[c(matrix(1:(n-2), byrow = TRUE, ncol = 22)[, 1], n), ]
        Table = Table[-(1:17),]
        Table = t(Table)
    } else if (size == "mini") {
        p = c(
            "0.01", "0.0251189", "0.05", "0.1",
            "0.9", "0.95", "0.9751687", "0.99")
        Table = Table[, p]
        Table = rbind(Table, qchisq(as.numeric(p), df = 2))
        Table = round(Table, 3)
        colnames(Table) = substr(colnames(Table), 1, 5)
        nRows = dim(Table)[1]
        rownames(Table)[nRows] = "Inf"
    }
    ans = list(
        x = as.numeric(rownames(Table)),
        y = as.numeric(colnames(Table)),
        z = Table)

    # Add Control:
    attr(ans, "control") <-
        c(table = "jb", type = type, size = size)

    # Return Value:
    class(ans) = "gridData"
    ans
}


# ------------------------------------------------------------------------------


.jbPlot <- 
function(type = c("LM", "ALM"))
{   
    # A function implemented by Diethelm Wuertz

    # Match Arguments:
    type = match.arg(type)

    # Load Table:
    Y = .jbTable(size = "small")
    X = cbind(expand.grid(x = Y$x, y = Y$y), z = as.vector(Y$z))
    x = X[, 1] # N
    y = X[, 3] # q-Stat
    z = X[, 2] # p-Value

    # Interpolate:
    ans = akimaInterp(log(x), log(y), z, gridPoints = 101, extrap = FALSE)
    persp(ans, theta = 40, xlab = "log(N)", ylab = "log(q)", zlab = "p",
        main = paste(type, "Jarque Bera"))

    # Return Value:
    invisible(NULL)
}



# ------------------------------------------------------------------------------


.pjb <- 
function(q, N = Inf, type = c("LM", "ALM"))
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Returns probabilities for the JB Test given quantiles

    # Arguments:
    #   q = quantiles (0.020 < q < 12)
    #   N = sample sizes

    # Details:
    #   Uses "small" table for interpolation.

    # Example:
    #   q=c(0.03,1,5,9); for (N in c(4,1000,50000,Inf)) print(.pjb(q, N))

    # FUNCTION:

    # Match Arguments:
    type = match.arg(type)

    # Check N:
    stopifnot(length(N) == 1)

    # Grid Points:
    Y = .jbTable(type = type, size = "small")
    if (N < 10000) {
        # Interpolate on log-Scale:
        X = cbind(expand.grid(x = Y$x, y = Y$y), z = as.vector(Y$z))
        x = log(X[, 1]) # N
        y = log(X[, 3]) # q-Stat
        z = X[, 2] # p-Value
        # Out Points:
        xo = rep(log(N), times = length(q))
        yo = log(q)
        # Interpolate:
        p = linearInterpp(x, y, z, xo = xo, yo = yo)[, 3]
    } else if (N > 10000) {
        p = pchisq(as.numeric(q), df = 2)
    }

    # Result:
    if (N < 5) controlN = "< 5"
    if (N >= 5 & N <= 10000) controlN = as.character(N)
    if (N > 10000) controlN = "> 10000"
    if (!is.finite(N)) controlN = "Inf"
    names(p) = as.character(q)
    attr(p, "control") <- c(N = controlN)

    # Return Value:
    p
}


# ------------------------------------------------------------------------------


.qjb <- 
function(p, N = Inf, type = c("LM", "ALM"))
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Returns quantiles for the ADF Test given probabilities

    # Arguments:
    #   p = probabilities (0.01 < p < 0.99)
    #   N = sample sizes

    # Details:
    #   Uses "small" table for interpolation.

    # Example:
    #   p=(1:9)/10; for (N in c(4,1000,50000,Inf)) print(.qjb(p, N))

    # FUNCTION:

    # Match Arguments:
    type = match.arg(type)
    statistic = match.arg(statistic)

    # Check N:
    stopifnot(length(N) == 1)

    # Grid Points:
    Y = .jbTable(type = type, size = "small")
    if (N < 10000) {
        # Interpolate on log-Scale:
        X = cbind(expand.grid(x = Y$x, y = Y$y), z = as.vector(Y$z))
        x = log(X[, 1]) # N
        y = X[, 2] # p-Value
        z = X[, 3] # q-Stat
        # Out Points:
        xo = rep(log(N), times = length(p))
        yo = p
        # Interpolate - Scaling Required:
        xStd = sqrt(var(x))
        q = linearInterpp(x/xStd, y, z, xo = xo/xStd, yo = yo)[, 3]
    } else {
        # Asymptotic Limit:
        q = qchisq(as.numeric(p), df = 2)
    }

    # Result:
    if (N < 5) controlN = "< 5"
    if (N >= 5 & N <= 10000) controlN = as.character(N)
    if (N > 10000) controlN = "> 10000"
    if (!is.finite(N)) controlN = "Inf"
    names(q) = as.character(p)
    attr(q, "control") <- c(N = controlN)

    # Return Value:
    q
}


################################################################################

