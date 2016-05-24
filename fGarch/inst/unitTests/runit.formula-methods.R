
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
#   1999 - 2008, Diethelm Wuertz, Rmetrics Foundation, GPL
#   Diethelm Wuertz <wuertz@itp.phys.ethz.ch>
#   info@rmetrics.org
#   www.rmetrics.org
# for the code accessed (or partly included) from other R-ports:
#   see R's copyright and license files
# for the code accessed (or partly included) from contributed R-ports
# and other sources
#   see Rmetrics's copyright file


################################################################################


test.formula.methods.univariate <-
    function()
{

    # Numeric Vector RVs:
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")

    # Simulate normal GARCH(1, 1) numeric Vector:
    spec = garchSpec()
    N = 250

    # Univariate Data Simulation:
    x.vec = 100*garchSim(spec, N)
    print(head(x.vec))
    x.tS = dummyDailySeries(matrix(x.vec), units = "GARCH11")
    print(head(x.tS))
#    x.zoo = zoo(as.vector(x.vec), order.by = as.Date(rownames(x.tS)))
#    print(head(x.zoo))
    x.ts = as.ts(x.vec)
    print(head(x.ts))

    # Univariate Modeling:

    # A numeric Vector:
    fit = garchFit(~ garch(1,1), data = x.vec, trace = FALSE)
    print(formula(fit))
    fit = garchFit(x.vec ~ garch(1,1), data = x.vec, trace = FALSE)
    print(formula(fit))

    # An univariate timeSeries object with dummy dates:
    fit = garchFit(~ garch(1,1), data = x.tS, trace = FALSE)
    print(formula(fit))
    fit = garchFit(x.tS ~ garch(1,1), data = x.tS, trace = FALSE)
    print(formula(fit))

###     # An univariate zoo object with dummy dates:
###     fit = garchFit(~ garch(1,1), data = x.zoo, trace = FALSE)
###     print(formula(fit))
###     fit = garchFit(x.zoo ~ garch(1,1), data = x.zoo, trace = FALSE)
###     print(formula(fit))

    # An univariate "ts" object:
    fit = garchFit(~ garch(1,1), data = x.ts, trace = FALSE)
    print(formula(fit))
    fit = garchFit(x.ts ~ garch(1,1), data = x.ts, trace = FALSE)
    print(formula(fit))

    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.formula.methods.multivariate <-
    function()
{
    # Numeric Vector RVs:
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")

    # Simulate normal GARCH(1, 1) numeric Vector:
    spec = garchSpec()
    N = 250

    # Univariate Data Simulation:
    x.vec = 100*garchSim(spec, N)
    print(head(x.vec))
    x.tS = dummyDailySeries(matrix(x.vec), units = "GARCH11")
    print(head(x.tS))

    # Multivariate Data Simulation:
    X.mat = cbind(GARCH11 = x.vec, R = rnorm(N))
    colnames(X.mat) <- c("GARCH11", "R")
    print(head(X.mat))
    X.tS = dummyDailySeries(X.mat, units = c("GARCH11", "R"))
    print(head(X.tS))
#    X.zoo = zoo(X.mat, order.by = as.Date(rownames(x.tS)))
#    print(head(X.zoo))
    X.mts = as.ts(X.mat)
    print(head(X.mts)) # head doesn't wor for mts !!!

    # Multivariate Modeling:

    # A numeric matrix:
    fit = garchFit(GARCH11 ~ garch(1,1), data = X.mat, trace = FALSE)
    print(formula(fit))
    fit = garchFit(100*GARCH11 ~ garch(1,1), data = X.mat, trace = FALSE)
    print(formula(fit))

    # A multivariate timeSeries object with dummy dates:
    fit = garchFit(GARCH11 ~ garch(1,1), data = X.tS, trace = FALSE)
    print(formula(fit))
    fit = garchFit(100*GARCH11 ~ garch(1,1), data = X.tS, trace = FALSE)
    print(formula(fit))

###     # A multivariate zoo object without column names:
###     fit = garchFit(GARCH11 ~ garch(1,1), data = X.zoo, trace = FALSE)
###     print(formula(fit))
###     fit = garchFit(100*GARCH11 + R/100 ~ garch(1,1), data = X.zoo, trace = FALSE)
###     print(formula(fit))

    # A multivariate "mts" object without column names:
    fit = garchFit(GARCH11 ~ garch(1,1), data = X.mts, trace = FALSE)
    print(formula(fit))
    fit = garchFit(100*GARCH11 + R/100 ~ garch(1,1), data = X.mts, trace = FALSE)
    print(formula(fit))

    # Return Value:
    return()
}

# ------------------------------------------------------------------------------


test.formula.methods.spread <-
    function()
{
    # MODELING THE PERCENTUAL SPI/SBI SPREAD FROM LPP BENCHMARK:

    # Series:
    X.tS = as.timeSeries(data(LPP2005REC))
    print(head(X.tS))
    X.mat = as.matrix(X.tS)
    print(head(X.mat))
#    X.zoo = zoo(X.mat, order.by = as.Date(rownames(X.tS)))
#    print(head(X.zoo))
    X.mts = ts(X.mat)
    print(head(X.mts)) # head does not work for ts objects!

    # Fit:
    fit = garchFit(100*(SPI - SBI) ~ garch(1,1), data = X.tS, trace = FALSE)
    print(formula(fit))
    ## fit = garchFit(100*(SPI - SBI) ~ garch(1,1), data = X.mat, trace = FALSE)
    ## print(formula(fit))
    ## fit = garchFit(100*(SPI - SBI) ~ garch(1,1), data = X.zoo, trace = FALSE)
    ## print(formula(fit))
    ## fit = garchFit(100*(SPI - SBI) ~ garch(1,1), data = X.mts, trace = FALSE)
    ## print(formula(fit))

    # MODELING HIGH/LOW SPREADS FROM MSFT PRICE SERIES:

    # Series:
    X.tS = MSFT

    # Fit:
    fit = garchFit(Open ~ garch(1,1), data = returns(X.tS), trace = FALSE)
    print(formula(fit))
    fit = garchFit(100*(High-Low) ~ garch(1,1), data = returns(X.tS),
        trace = FALSE)
    print(formula(fit))

    # Return Value:
    return()
}


################################################################################

