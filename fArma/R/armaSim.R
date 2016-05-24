
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
# FUNCTION:               SIMULATION:
#  armaSim                 Simulates an ARIMA time series process
################################################################################


armaSim <-
    function(model = list(ar = c(0.5, -0.5), d = 0, ma = 0.1), n = 100,
             innov = NULL, n.start = 100, start.innov = NULL,
             rand.gen = rnorm, rseed = NULL, addControl = FALSE, ...)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Simulates an ARIMA Time Series Process

    # Details:
    #   Splus-Like argument list ...
    #   Rmetrics Notation:
    #     armaSim(model = list(ar = c(0.5, -0.5), d = 0, ma = 0.1), n = 100,
    #       innov = NULL, n.start = 100, start.innov = NULL,
    #       rand.gen = rnorm, rseed = NULL, ...)
    # SPlus Notation:
    #     arima.sim (model, n = 100,
    #       innov = rand.gen(n, ...), n.start = 100, start.innov = NULL,
    #       rand.gen = rnorm, xreg = NULL, reg.coef = NULL, ...)

    # Example:
    #   armaSim(model = list(ar = c(0.5, -0.5), d = 0, ma = 0.1))
    #   armaSim(model = list(ar = c(0.5, -0.5), d = 0.2, ma = 0.1))
    #   armaSim(model = list(ar = 0, d = 0.2, ma = 0))
    #   armaSim(model = list(d = 0.2))

    # FUNCTION:

    # Checks:
    if (!is.list(model))
        stop("model must be a list")

    # Simulate:
    if (!is.null(rseed))
        set.seed(rseed)
    if (is.null(innov))
        innov = rand.gen(n, ...)
    n = length(innov)
    if (is.null(start.innov))
        start.innov = rand.gen(n, ...)
    n.start = length(start.innov)

    # AR PART:
    p = length(model$ar)
    if (p == 1 && model$ar == 0)
        p = 0
    if (p) {
        minroots = min(Mod(polyroot(c(1, -model$ar))))
        if (minroots <= 1) warning(" AR part of model is not stationary")
    }

    # MA PART:
    q = length(model$ma)
    if (q == 1 && model$ma == 0)
        q = 0
    if (n.start < p + q)
        stop("burn-in must be as long as ar + ma")

    # DIFFERENCING:
    ## if (model$d < 0) stop("d must be positive ")
    dd = length(model$d)
    if (dd) {
        # ARFIMA|FRACDIFF if "dd" is a non-integer value:
        d = model$d
        if (d != round(d) ) {
            TSMODEL = "ARFIMA"
        } else {
            TSMODEL = "ARIMA"
        }
    } else {
        d = 0
        TSMODEL = "ARIMA"
    }

    # ARMA:
    if (TSMODEL == "ARIMA") {
        x = ts(c(start.innov, innov), start = 1 - n.start)
        if (length(model$ma)) {
            x <- filter(x, c(1, model$ma), sides = 1)
            x[seq_along(model$ma)] <- 0
        }
        if (length(model$ar)) x = filter(x, model$ar, method = "recursive")
        x = x[-(1:n.start)]
        if (d > 0) x = diffinv(x, differences = d)
    }

    # ARFIMA [FRACDIFF]:
    if (TSMODEL == "ARFIMA") {
        if (p == 0) model$ar = 0
        if (q == 0) model$ma = 0
        mu = 0
        # Use Fortran Routine from R's contributed fracdiff package:
        # This is a BUILTIN function ...
        if (!is.null(rseed)) set.seed(rseed)
        eps = rnorm(n + q)
        x = .Fortran("fdsim", as.integer(n), as.integer(p), as.integer(q),
            as.double(model$ar), as.double(model$ma), as.double(model$d),
            as.double(mu), as.double(eps), x = double(n + q),
            as.double(.Machine$double.xmin), as.double(.Machine$double.xmax),
            as.double(.Machine$double.neg.eps), as.double(.Machine$double.eps),
            PACKAGE = "fArma")$x[1:n]
    }

    # Time Series:
    from <-
        timeDate(format(Sys.time(), format = "%Y-%m-%d")) - NROW(x)*24*3600
    charvec  <- timeSequence(from = from, length.out = NROW(x))
    ans = timeSeries(data = matrix(x, ncol = 1), charvec = charvec, ...)

    # Add Control:
    if (addControl) {
        control = c(ar = model$ar, d = model$d, ma = model$ma)
        Names = names(control)
        control = as.character(c(control, substitute(rand.gen)))
        names(control) = c(Names, "rand.gen")
        if (!is.null(rseed)) control = c(control, rseed = rseed)
        attr(ans, "control") = control
    }

    # Return Value:
    ans
}


################################################################################

