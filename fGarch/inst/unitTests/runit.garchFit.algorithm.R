
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


test.garchFit.nlminb <-
    function()
{
    # RVs:
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")

    # Load Data:
    x = garchSim(n = 100)

    # Algorithms:
    #   "nlminb", "sqp", "lbfgsb", "nlminb+nm", "lbfgsb+nm"

    # nlminb:
    fit = garchFit( ~ garch(1,1), data = x,
        algorithm = "nlminb", trace = FALSE)
    print(coef(fit))

    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.garchFit.sqp <-
    function()
{
    # RVs:
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")

    # Load Data:
    x = garchSim(n = 100)

    # Algorithms:
    #   "nlminb", "sqp", "lbfgsb", "nlminb+nm", "lbfgsb+nm"

###     # sqp:
###     fit = garchFit( ~ garch(1,1), data = x,
###         algorithm = "sqp", trace = FALSE)
###     print(coef(fit))


    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.garchFit.lbfgsb <-
    function()
{
    # RVs:
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")

    # Load Data:
    x = garchSim(n = 100)

    # Algorithms:
    #   "nlminb", "sqp", "lbfgsb", "nlminb+nm", "lbfgsb+nm"

    # lbfgsb:
    fit = garchFit( ~ garch(1,1), data = x,
        algorithm = "lbfgsb", trace = FALSE)
    coef(fit)

    # nlminb+nm:
    fit = garchFit( ~ garch(1,1), data = x,
        algorithm = "nlminb+nm", trace = FALSE)
    coef(fit)

    # lbfgsb+nm:
    fit = garchFit( ~ garch(1,1), data = x,
        algorithm = "lbfgsb+nm", trace = FALSE)
    coef(fit)

    # Return Value:
    return()
}

# ------------------------------------------------------------------------------


test.garchFit.nlmin.nm <-
    function()
{
    # RVs:
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")

    # Load Data:
    x = garchSim(n = 100)

    # Algorithms:
    #   "nlminb", "sqp", "lbfgsb", "nlminb+nm", "lbfgsb+nm"

    # nlminb+nm:
    fit = garchFit( ~ garch(1,1), data = x,
        algorithm = "nlminb+nm", trace = FALSE)
    coef(fit)

    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.garchFit.lbfgsb.nm <-
    function()
{
    # RVs:
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")

    # Load Data:
    x = garchSim(n = 100)

    # Algorithms:
    #   "nlminb", "sqp", "lbfgsb", "nlminb+nm", "lbfgsb+nm"

    # lbfgsb+nm:
    fit = garchFit( ~ garch(1,1), data = x,
        algorithm = "lbfgsb+nm", trace = FALSE)
    coef(fit)

    # Return Value:
    return()
}

################################################################################

