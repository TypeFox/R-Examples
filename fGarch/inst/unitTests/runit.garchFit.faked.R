
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
# FUNCTION:               PARAMETER ESTIMATION:
#  'fGARCH'                S4: fGARCH Class representation
#  garchFit                Fits GARCH and APARCH processes
################################################################################


test.garchFit.faked <-
function()
{
    # Numeric Vector RVs:
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")

    # Simulate normal GARCH(1, 1) numeric Vector:
    model = list(omega = 1e-06, alpha = 0.1, beta = 0.8)
    spec = garchSpec(model)
    print(spec)
    N = 250

    # UNIVARIATE:

    x.vec = as.vector(100*garchSim(spec, N))
    print(head(x.vec))
    x.tS = dummyDailySeries(matrix(x.vec), units = "GARCH11")
    print(head(x.tS))
#    x.zoo = zoo(as.vector(x.vec), order.by = as.Date(rownames(x.tS)))
#    print(head(x.zoo))
    x.ts = ts(x.vec)
    print(head(x.ts))

    # FIT FROM UNIVARIATE DATA SERIES:

    fit = garchFit( ~ garch(1,1), data = x.tS, trace = FALSE)
    print(fit)
    formula(fit)

    fit = garchFit( ~ garch(1,1), data = as.vector(x.tS),
        trace = FALSE)
    print(fit)
    formula(fit)

    a = 2
    b = 2
    fit = garchFit( ~ garch(1,1), data = a*as.vector(0+b*x.tS),
        trace = FALSE)
    print(fit)
    formula(fit)

    # WHAT HAPPENS WHEN WE (MIS)SPECIFY LHS ?
    #   ... lhs will be ignored for the univariate case:
    fit = garchFit(any ~ garch(1,1), data = x.vec, trace = FALSE)
    print(fit)
    formula(fit)

    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.garchFit.mult.faked <-
function()
{
    # Numeric Vector RVs:
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")

    # Simulate normal GARCH(1, 1) numeric Vector:
    model = list(omega = 1e-06, alpha = 0.1, beta = 0.8)
    spec = garchSpec(model)
    print(spec)
    N = 250

    # UNIVARIATE:

    x.vec = as.vector(100*garchSim(spec, N))
    print(head(x.vec))
    x.tS = dummyDailySeries(matrix(x.vec), units = "GARCH11")
    print(head(x.tS))
#    x.zoo = zoo(as.vector(x.vec), order.by = as.Date(rownames(x.tS)))
#    print(head(x.zoo))
    x.ts = ts(x.vec)
    print(head(x.ts))

    # MULTIVARIATE:

    X.mat = cbind(GARCH11 = x.vec, R = rnorm(N)/1000)
    print(head(X.mat))
    X.tS = dummyDailySeries(X.mat, units = c("GARCH11", "R"))
    print(head(X.tS))
#    X.zoo = zoo(X.mat, order.by = as.Date(rownames(x.tS)))
#    print(head(X.zoo))
    X.mts = ts(X.mat)
    print(head(X.mts))

    # FIT FROM MULTIVARIATE DATA SET:

    fit = garchFit(GARCH11 ~ garch(1,1), data = X.tS, trace = FALSE)
    print(fit)
    formula(fit)

    fit = garchFit(GARCH11 ~ garch(1,1), data = as.matrix(X.tS),
        trace = FALSE)
    print(fit)
    formula(fit)

    a = 2
    b = 2
    fit = garchFit(GARCH11 ~ garch(1,1), data = a*as.matrix(0+b*X.tS),
        trace = FALSE)
    print(fit)
    formula(fit)

    a = 2
    b = 2
    fit = garchFit(GARCH11 ~ garch(1,1), data = a*as.matrix(0+b*X.tS),
        trace = FALSE)
    print(fit)
    formula(fit)

    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.garchFit.mult.lhs.faked <-
function()
{
    # Numeric Vector RVs:
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")

    # Simulate normal GARCH(1, 1) numeric Vector:
    model = list(omega = 1e-06, alpha = 0.1, beta = 0.8)
    spec = garchSpec(model)
    print(spec)
    N = 250

    # UNIVARIATE:

    x.vec = as.vector(100*garchSim(spec, N))
    print(head(x.vec))
    x.tS = dummyDailySeries(matrix(x.vec), units = "GARCH11")
    print(head(x.tS))
#    x.zoo = zoo(as.vector(x.vec), order.by = as.Date(rownames(x.tS)))
#    print(head(x.zoo))
    x.ts = ts(x.vec)
    print(head(x.ts))

    # MULTIVARIATE:

    X.mat = cbind(GARCH11 = x.vec, R = rnorm(N)/1000)
    print(head(X.mat))
    X.tS = dummyDailySeries(X.mat, units = c("GARCH11", "R"))
    print(head(X.tS))
#    X.zoo = zoo(X.mat, order.by = as.Date(rownames(x.tS)))
#    print(head(X.zoo))
    X.mts = ts(X.mat)
    print(head(X.mts))

    # LEFT HAND SIDE FORMULA FAKED FIT:

    fit = garchFit(GARCH11 + R ~ garch(1,1), data = X.tS, trace = FALSE)
    print(fit)
    formula(fit)
###     head(fit@data$data)
###     head(fit@data$Data)
###     head(rowSums(fit@data$Data))

    # LEFT HAND SIDE FORMULA FAKED AND DATA FAKED FIT:

    fit = garchFit(GARCH11 + R ~ garch(1,1), data = as.matrix(X.tS),
        trace = FALSE)
    print(fit)
    formula(fit)
###     head(fit@data$data)
###     head(fit@data$Data)
###     head(rowSums(fit@data$Data))

    # Return Value:
    return()
}


################################################################################

