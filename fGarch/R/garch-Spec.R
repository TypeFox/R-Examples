
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
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
# FUNCTION:               SPECIFICATION:
#  garchSpec               Creates a 'garchSpec' object from scratch
###############################################################################


garchSpec <-
    function (model = list(), presample = NULL,
    cond.dist = c("norm", "ged", "std", "snorm", "sged", "sstd"),
    rseed = NULL)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Creates a "garchSpec" object from scratch.

    # Arguments:
    #   model - a list with the model parameters as entries
    #     omega - the variance value for GARCH/APARCH
    #       specification,
    #     alpha - a vector of autoregressive coefficients
    #       of length p for the GARCH/APARCH specification,
    #     gamma - a vector of leverage coefficients of
    #       length p for the APARCH specification,
    #     beta - a vector of moving average coefficients of
    #       length q for the GARCH/APARCH specification,
    #     mu - the mean value for ARMA specification,
    #     ar - a vector of autoregressive coefficients of
    #       length m for the ARMA specification,
    #     ma - a vector of moving average coefficients of
    #       length n for the ARMA specification,
    #     delta - the exponent value used in the variance equation.
    #     skew - a numeric value listing the distributional
    #        skewness parameter.
    #     shape - a numeric value listing the distributional
    #        shape parameter.
    #   presample - either a multivariate "timeSeries", a
    #       multivariate "ts", a "data.frame" object or a numeric
    #       "matrix" with 3 columns and at least max(m,n,p,q)
    #       rows. The first culumn are the innovations, the second
    #       the conditional variances, and the last the time series.
    #   condd.dist - a character string naming the distribution
    #       function.
    #   rseed - optional random seed.

    # Slots:
    #   call - the function call.
    #   formula - a formula object describing the model, e.g.
    #       ARMA(m,n) + GARCH(p,q). ARMA can be missing or
    #       specified as AR(m) or MA(n) in the case of pure
    #       autoregressive or moving average models. GARCH may
    #       alternatively specified as ARCH(p) or APARCH(p,q).
    #       If formula is set to "NA", the formula is constructed
    #       from the "model" list.
    #   model - as declared in the input.

    # FUNCTION:

    # Match Arguments:
    cond.dist = match.arg(cond.dist)

    # Skewness Parameter Settings:
    skew = list(
        "norm" = NULL,
        "ged" = NULL,
        "std" = NULL,
        "snorm" = 0.9,
        "sged" = 0.9,
        "sstd" = 0.9)

    # Shape Parameter Settings:
    shape = list(
        "norm" = NULL,
        "ged" = 2,
        "std" = 4,
        "snorm" = NULL,
        "sged" = 2,
        "sstd" = 4)

    # Default Model:
    control = list(
        omega = 1.0e-6,
        alpha = 0.1,
        gamma = NULL,
        beta = 0.8,
        mu = NULL,
        ar = NULL,
        ma = NULL,
        delta = 2,
        skew = skew[[cond.dist]],
        shape = shape[[cond.dist]]
        )

    # Update Control:
    control[names(model)] <- model
    model <- control

    # check if alpha and beta are well defined
    if (sum(c(model$alpha, model$beta))>1)
        warnings("sum(alpha)+sum(beta)>1")

    # Model Orders:
    order.ar = length(model$ar)
    order.ma = length(model$ma)
    order.alpha = length(model$alpha)
    if (sum(model$beta) == 0) {
        order.beta = 0
    } else {
        order.beta = length(model$beta)
    }

    # Compose Mean Formula Object:
    if (order.ar == 0 && order.ma == 0) {
        formula.mean = ""
    }
    if (order.ar > 0 && order.ma == 0) {
        formula.mean = paste ("ar(", as.character(order.ar), ")", sep = "")
    }
    if (order.ar == 0 && order.ma > 0) {
        formula.mean = paste ("ma(", as.character(order.ma), ")",  sep = "")
    }
    if (order.ar > 0 && order.ma > 0) {
        formula.mean = paste ("arma(", as.character(order.ar), ", ",
            as.character(order.ma), ")", sep = "")
    }

    # Compose Variance Formula Object:
    formula.var = "garch"
    if (order.beta == 0) formula.var = "arch"
    if (!is.null(model$gamma) != 0) formula.var = "aparch"
    if (model$delta != 2) formula.var = "aparch"

    if (order.beta == 0) {
        formula.var = paste(formula.var, "(", as.character(order.alpha), ")",
            sep = "")
    } else {
        formula.var = paste(formula.var, "(", as.character(order.alpha),
            ", ", as.character(order.beta), ")", sep = "")
    }

    # Compose Mean-Variance Formula Object:
    if (formula.mean == "") {
        formula = as.formula(paste("~", formula.var))
    } else {
        formula = as.formula(paste("~", formula.mean, "+", formula.var))
    }

    # Add NULL default entries:
    if (is.null(model$mu)) model$mu = 0
    if (is.null(model$ar)) model$ar = 0
    if (is.null(model$ma)) model$ma = 0
    if (is.null(model$gamma)) model$gamma = rep(0, times = order.alpha)
    # print(unlist(model))

    # Seed:
    if (is.null(rseed)) {
        rseed = 0
    } else {
        set.seed(rseed)
    }

    # Define Missing Presample:
    order.max = max(order.ar, order.ma, order.alpha, order.beta)
    iterate = TRUE
    if (!is.matrix(presample)) {
        if (is.null(presample)) {
            iterate = FALSE
            n.start = order.max
        } else {
            n.start = presample
        }
        z = rnorm(n = n.start)
        # GARCH(p, q):
        h = rep(model$omega/(1-sum(model$alpha)-sum(model$beta)),
            times = n.start)
        y = rep(model$mu/(1-sum(model$ar)), times = n.start)
        # APARCH(p,q):
        #  ... we initialize all models with norm-GARCH(p,q) processes
    } else {
        z = presample[, 1]
        h = presample[, 2]
        y = presample[, 3]
    }
    presample = cbind(z, h, y)

    # Presample Iteration:
    if (iterate) {
        n.iterate = length(z) - order.max
        deltainv = 1/model$delta
        for (i in n.iterate:1) {
            h[i] = model$omega +
                sum(model$alpha*(abs(abs(y[i+(1:order.alpha)]) -
                    model$gamma*y[i+(1:order.alpha)])^model$delta)) +
                sum(model$beta*h[i+(1:order.beta)])
            y[i] = model$mu  +
                sum(model$ar*y[i+(1:order.ar)]) +
                sum(model$ma*(h[i+(1:order.ma)]**deltainv)) +
                h[i]^deltainv * z[i]
        }
    }

    # Result:
    new("fGARCHSPEC",
        call = match.call(),
        formula = formula,
        model = list(omega = model$omega, alpha = model$alpha,
            gamma = model$gamma, beta = model$beta, mu = model$mu,
            ar = model$ar, ma = model$ma, delta = model$delta,
            skew = model$skew, shape = model$shape),
        presample = as.matrix(presample),
        distribution = as.character(cond.dist),
        rseed = as.numeric(rseed)
    )
}


################################################################################

