
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


# garchFit(
    #   formula,
    #   data,
    #   init.rec = c("mci", "uev"),
    #   delta = 2,
    #   skew = 1,
    #   shape = 4,
    #   cond.dist = c("dnorm", "dsnorm", "dged", "dsged", "dstd", "dsstd"),
    #   include.mean = TRUE,
    #   include.delta = NULL,
    #   include.skew = NULL,
    #   include.shape = NULL,
    #   leverage = NULL,
    #   trace = TRUE,
    #   algorithm = c("sqp", "nlminb", "lbfgsb", "nlminb+nm", "lbfgsb+nm"),
    #   control = list(),
    #   title = NULL,
    #   description = NULL,
    #   ...)


# ------------------------------------------------------------------------------


test.garchInputSeries <-
function()
{
    # Numeric Vector RVs:
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")

    # Simulate normal GARCH(1, 1) numeric Vector:
    model = list(omega = 1e-06, alpha = 0.1, beta = 0.8)
    spec = garchSpec(model)
    print(spec)
    N = 10


    # UNIVARIATE:

    # A numeric Vector:
    x.vec = as.vector(100*garchSim(spec, N))
    print(head(x.vec))
    x.tS = dummyDailySeries(matrix(x.vec), units = "GARCH11")
    print(head(x.tS))
#   x.zoo = zoo(as.vector(x.vec), order.by = as.Date(rownames(x.tS)))
#    print(head(x.zoo))
    x.ts = ts(x.vec)
    print(head(x.ts))

    # MULTIVARIATE:

    # A numeric matrix:
    X.mat = cbind(GARCH11 = x.vec, R = rnorm(N))
    print(head(X.mat))
    X.tS = dummyDailySeries(X.mat, units = c("GARCH11", "R"))
    print(head(X.tS))
#    X.zoo = zoo(X.mat, order.by = as.Date(rownames(x.tS)))
#    print(head(X.zoo))
    X.mts = ts(X.mat)
    print(head(X.mts))

    # Return Value:
    return()
}


################################################################################

