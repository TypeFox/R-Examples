
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


test.garchFit.garch11 <-
    function()
{
    # Use Simulated Series - an Object of class 'ts' ...

    # RVs:
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")

    # Normal GARCH(1, 1)
    x = garchSim(n = 250)

    # Fit:
    fit = garchFit( ~ garch(1,1), data = x, trace = FALSE)
    print(coef(fit))

    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.garchFit.garch21 <-
    function()
{
    # Use Simulated Series - an Object of class 'ts' ...

    # RVs:
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")

    # Normal-GARCH(2, 1)
    model = list(omega = 1e-06, alpha = c(0.1, 0.2), beta = 0.6)
    spec = garchSpec(model)
    x = garchSim(spec = spec, n = 250)

    # Fit
    fit = garchFit( ~ garch(2,1), data = x, trace = FALSE)
    print(coef(fit))

    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.garchFit.ar1garch11 <-
    function()
{
    # Use Simulated Series - an Object of class 'ts' ...

    # RVs:
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")

    # Normal AR(1)-GARCH(1,1):
    model = list(omega = 1e-06, ar = -0.1, alpha = c(0.1, 0.2), beta = 0.6)
    spec = garchSpec(model)
    x = garchSim(spec = spec, n = 250)

    # Fit:
    fit = garchFit(~ arma(1,0) + garch(1,1), data = x, trace = FALSE)
    print(coef(fit))

    # Return Value:
    return()
}


################################################################################

