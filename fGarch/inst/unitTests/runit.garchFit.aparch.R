
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


test.aparch11 <-
    function()
{
    # Use Simulated Series - an Object of class 'ts' ...

    # RVs:
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")

    # Leveraged normal APARCH(1, 1)
    model = list(omega = 1e-06, alpha = 0.1, beta = 0.8, gamma = 0.1)
    spec = garchSpec(model)
    x = garchSim(spec, n = 250)

    # Fit:
    fit = garchFit(x ~ garch(1,1), data = x, leverage = TRUE, trace = FALSE)
    print(coef(fit))

    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.aparch11delta <-
    function()
{
    # RVs:
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")

    # Leveraged normal APARCH(1, 1) delta = 1.5
    model = list(omega = 1e-05, alpha = 0.1, beta = 0.8, gamma = 0.1,
        delta = 1.5)
    spec = garchSpec(model)
    x = garchSim(spec, n = 250)
    print(var(x))

    # Fit:
    fit = garchFit(x ~ garch(1,1), data = x, leverage = TRUE,
        include.delta = TRUE, delta = 2, trace = FALSE)
    print(coef(fit))

    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.ar1aparch21 <-
    function()
{
    # RVs:
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")

    # Leveraged normal AR(1)-APARCH(2, 1) delta = 1
    model = list(omega = 1e-06, ar = -0.1, alpha = c(0.1, 0.2), beta = 0.6,
        delta = 1)
    spec = garchSpec(model)
    x = garchSim(spec, n = 250)

    # Taylor Plot:
    taylor = teffectPlot(x)
    init.delta = mean(taylor$maxDelta)
    init.delta

    ## fit = garchFit(x ~ ar(1) + garch(2,1), data = x, include.delta = TRUE,
    ##     delta = init.delta, trace = FALSE)
    ## print(coef(fit))

    ## Error in solve.default(fit$hessian) :
    ##  Lapack routine dgesv: system is exactly singular            ## CHECK !!!

    # Return Value:
    return()
}


################################################################################

