
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
#   1999 - 2007, Diethelm Wuertz, GPL
#   Diethelm Wuertz <wuertz@itp.phys.ethz.ch>
#   info@rmetrics.org
#   www.rmetrics.org
# for the code accessed (or partly included) from other R-ports:
#   see R's copyright and license files
# for the code accessed (or partly included) from contributed R-ports
# and other sources
#   see Rmetrics's copyright file


################################################################################
# FUNCTION:             GEV SIMULATION:
#  gevSim                Simulates a GEV distributed process
#  gumbelSim             Simulates a Gumbel distributed process
# FUNCTION:             GEV PARAMETER ESTIMATION:
#  'fGEVFIT'             S4 class representation
#  gevFit                Fits Parameters of GEV distribution
#  gumbelFit             Fits Parameters of Gumbel distribution
# METHODS:              PRINT, PLOT, AND SUMMARY:
#  show.fGEVFIT          S4 Show method for object of class "fGEVFIT"
#  plot.fGEVFIT          S3 Plot method for object of class "fGEVFIT"
#  summary.fGEVFIT       S3 Summary Method for object of class "fGEVFIT"
################################################################################


test.gevSim =
function()
{
    # gevSim(model = list(xi=-0.25, mu=0, beta=1), n = 1000, seed = NULL)

    # RVs:
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")

    # Artificial Data Set:
    model = list(xi = -0.25, mu = 0, beta = 1)
    x.ts = gevSim(model, n = 50, seed = 4711)
    class(x.ts)
    print(x.ts)

    # Create a daily timeSeries object with dummy dates:
    as.timeSeries(x.ts)

    # Create a daily timeSeries object starting 2007-01-01
    Calendar = timeSequence(from = "2007-01-01", length.out = length(x.ts))
    x.tS = timeSeries(data = x.ts, charvec = Calendar, units = "x")
    print(x.tS)

    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.gumbelSim =
function()
{
    # gumbelSim(model = list(mu=0, beta=1), n = 1000, seed = NULL)

    # RVs:
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")

    # Artificial Data Set:
    model = list(mu = 0, beta = 1)
    x.ts = gumbelSim(model, n = 50, seed = 4711)
    class(x.ts)
    print(x.ts)

    # Create a daily timeSeries object with dummy dates:
    x.tS = as.timeSeries(x.ts)
    print(x.tS)

    # Create a daily timeSeries object starting 2007-01-01
    Calendar = timeSequence(from = "2007-01-01", length.out = length(x.ts))
    x.tS = timeSeries(data = x.ts, charvec = Calendar, units = "x")
    print(x.tS)

    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.numericVectorBlocks =
function()
{
    # RVs:
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")

    # Check numeric vector as input:
    X = rt(5000, df = 4)
    x.vec = blockMaxima(X, 20)
    class(x.vec)
    head(x.vec)

    # Internal Fit - GEV PWM:
    fit = .gevpwmFit(x.vec)
    fit
    fit$par.ests

    # Internal Fit - GEV MLE:
    fit = .gevmleFit(x.vec)
    fit
    fit$par.ests

    # Internal Fit - Gumbel PWM:
    fit = .gumpwmFit(x.vec)
    fit
    fit$par.ests

    # Internal Fit - Gumbel MLE:
    fit = .gummleFit(x.vec)
    fit
    fit$par.ests

    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.timeSeriesBlocks =
function()
{
    # RVs:
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")

    # Create an artificial timeSeries with dummy positions:
    xx <- rt(5000, df = 4)
    charvec <- timeSequence("2013-01-01", length.out = NROW(xx))
    X = timeSeries(xx, charvec = charvec)

    # Compute Block Maxima:
    x.tS = blockMaxima(X, "monthly")
    class(x.tS)
    head(x.tS)

    # Convert to Vector:
    x.vec = as.vector(x.tS)

    # Internal Fit - GEV PWM:
    fit = .gevpwmFit(x.vec)
    fit
    fit$par.ests

    # Internal Fit - GEV MLE:
    fit = .gevmleFit(x.vec)
    fit
    fit$par.ests

    # Internal Fit - Gumbel PWM:
    fit = .gumpwmFit(x.vec)
    fit
    fit$par.ests

    # Internal Fit - Gumbel MLE:
    fit = .gummleFit(x.vec)
    fit
    fit$par.ests

    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.gevFit =
function()
{
    # gevFit(x, block = 1, type = c("mle", "pwm"),
    #   title = NULL, description = NULL, ...)

    # RVs:
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")

    # Simulate Series:
    model = list(xi = -0.25, mu = 0, beta = 1)
    x.ts = gevSim(model = model, n = 5000, seed = 4711)
    class(x.ts)

    # Check time series input:
    fit = gevFit(x.ts, block = 1, type = "pwm")
    class(fit)
    print(fit)
    fit = gevFit(x.ts, block = 1, type = "mle")
    class(fit)
    print(fit)

    # Check numeric vector input:
    fit = gevFit(as.vector(x.ts), block = 1, type = "pwm")
    class(fit)
    print(fit)
    fit = gevFit(as.vector(x.ts), block = 1, type = "mle")
    class(fit)
    print(fit)

    # Check timeSeries objerct input:
    fit = gevFit(as.timeSeries(x.ts), block = 1, type = "pwm")
    class(fit)
    print(fit)
    fit = gevFit(as.timeSeries(x.ts), block = 1, type = "mle")
    class(fit)
    print(fit)

    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.gumbelFit =
function()
{
    # gevFit(x, block = 1, type = c("mle", "pwm"),
    #   title = NULL, description = NULL, ...)

    # RVs:
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")

    # Simulate Series:
    model = list(mu = 0, beta = 1)
    x.ts = gumbelSim(model = model, n = 5000, seed = 4711)
    class(x.ts)

    # Check time series input:
    fit = gumbelFit(x.ts, block = 1, type = "pwm")
    class(fit)
    print(fit)
    fit = gumbelFit(x.ts, block = 1, type = "mle")
    class(fit)
    print(fit)

    # Check numeric vector input:
    fit = gumbelFit(as.vector(x.ts), block = 1, type = "pwm")
    class(fit)
    print(fit)
    fit = gumbelFit(as.vector(x.ts), block = 1, type = "mle")
    class(fit)
    print(fit)

    # Check timeSeries objerct input:
    fit = gumbelFit(as.timeSeries(x.ts), block = 1, type = "pwm")
    class(fit)
    print(fit)
    fit = gumbelFit(as.timeSeries(x.ts), block = 1, type = "mle")
    class(fit)
    print(fit)

    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.gevFitByBlocks <-
function()
{
    # gevFit(x, block = 1, type = c("mle", "pwm"),
    #   title = NULL, description = NULL, ...)

    # RVs:
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")

    # Simulate Series:
    model = list(xi = -0.25, mu = 0, beta = 1)
    x.ts = gevSim(model = model, n = 5000, seed = 4711)
    class(x.ts)
    x.vec = as.vector(x.ts)
    class(x.vec)
    charvec <- timeSequence("2013-01-01", length.out = NROW(x.vec))
    x.tS = timeSeries(x.vec, charvec)

    # ts as input & 20 Days Blocks:
    fit = gevFit(x.ts, block = 20, type = "pwm")
    fit
    fit = gevFit(x.ts, block = 20, type = "mle")
    fit

    # Numeric Vector as input & 20 Days Blocks:
    fit = gevFit(x.vec, block = 20, type = "pwm")
    fit
    fit = gevFit(x.vec, block = 20, type = "mle")
    fit

    # timeSeries o bject as input & Monthly Blocks:
    fit = gevFit(x.tS, block = "monthly", type = "pwm")
    fit
    fit = gevFit(x.tS, block = "quarterly", type = "mle")
    fit

    # timeSeries object as input & 20 Days Blocks:
    fit = gevFit(x.tS, block = 20, type = "pwm")
    fit
    fit = gevFit(x.tS, block = 20, type = "mle")
    fit

    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.plot =
function()
{
    # Load Data:
    x = as.timeSeries(data(danishClaims))

    # Parameter Estimation with Declustering:
    # gevFit(x, block = 1, type = c("mle", "pwm"),
    #   title = NULL, description = NULL, ...)
    fit = gevFit(x, block = "month")
    print(fit)

    # Plot:
    par(mfrow = c(2, 2), cex = 0.7)
    par(ask = FALSE)
    plot(fit, which = 1:4)

    # Try Interactive:
    # plot(fit)

    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.summary =
function()
{
    # Summary Report:
    # summary(object, doplot = TRUE, which = "all", ...)

    # Load Data:
    x = as.timeSeries(data(danishClaims))

    # Parameter Estimation with Declustering:
    fit = gevFit(x, block = "month")
    print(fit)

    # Summary:
    summary(fit, doplot = FALSE)

    # Return Value:
    return()
}


################################################################################
