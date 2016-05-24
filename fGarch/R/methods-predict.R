
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
# METHOD:                 PREDICTION:
#  predict.fGARCH          Forecasts from an object of class 'fGARCH'
################################################################################

setMethod(f = "predict", signature(object = "fGARCH"), definition =
          function(object, n.ahead = 10, trace = FALSE,
                   mse = c("cond","uncond"),
                   plot=FALSE, nx=NULL, crit_val=NULL, conf=NULL, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Prediction method for an object of class fGARCH

    # Arguments:
    #   object    an object of class fGARCH as returned by the
    #             function garchFit().
    #   n.ahead   number of steps to be forecasted, an integer
    #             value, by default 10)
    #   trace     should the prediction be traced? A logical value,
    #             by default FALSE)
    #    mse      should the mean squared errors be conditional or unconditional
    #    plot     should the predictions be plotted
    #    nx       The number of observations to be plotted with the predictions
    #             (If plot is TRUE, the default value of nx is the sample
    #             size times 0.25.)
    #    crit_va  If you want to set manually the critical values for
    #             the confidence intervals
    #    conf     The confidence level for computing the critical values
    #             of the confidence intervals

    # FUNCTION:

    mse <- match.arg(mse)

    # Retrieve "fit" from Parameter Estimation:
    fit = object@fit

    # Get ARMA(u,v)-GARCH(p,q) Order:
    u = fit$series$order[1]
    v = fit$series$order[2]
    p = fit$series$order[3]
    q = fit$series$order[4]
    max.order = max(u, v, p, q)

    # Get Start Conditions:
    h.start = fit$series$h.start
    llh.start = fit$series$llh.start
    index = fit$params$index
    params = fit$params$params
    par = fit$par
    Names = names(index)
    for (Name in Names) params[Name] = par[Name]
    Names = names(params)

    # Retrieve From Initialized Parameters:
    cond.dist = fit$params$cond.dist

    # Extract the Parameters by Name:
    leverage = fit$params$leverage
    mu = params["mu"]
    if (u > 0) {
        ar = params[substr(Names, 1, 2) == "ar"]
    } else {
        ar = c(ar1 = 0)
    }
    if (v > 0) {
        ma = params[substr(Names, 1, 2) == "ma"]
    } else {
        ma = c(ma1 = 0)
    }
    omega = params["omega"]
    if (p > 0) {
        alpha = params[substr(Names, 1, 5) == "alpha"]
    } else {
        alpha = c(alpha1 = 0)
    }
    if (p > 0 & leverage) {
        gamma = params[substr(Names, 1, 5) == "gamma"]
    } else {
        gamma = c(gamma1 = 0)
    }
    if (q > 0) {
        beta  = params[substr(Names, 1, 4) == "beta"]
    } else {
        beta = c(beta1 = 0)
    }
    delta = params["delta"]
    skew = params["skew"]
    shape = params["shape"]

    # Trace Parameters:
    if (trace) {
        cat("\nModel Parameters:\n")
        print(c(mu, ar, ma, omega, alpha, gamma, beta, delta, skew, shape))
    }

    # Retrieve Series Lengths:
    M = n.ahead
    N = length(object@data)

    # Get and Extend Series:
    x = c(object@data, rep(mu, M))
    h = c(object@h.t, rep(0, M))
    z = c(fit$series$z, rep(mu, M))

    # Forecast and Optionally Trace Variance Model:
    var.model = fit$series$model[2]
    # Forecast GARCH Variance:
    if (var.model == "garch") {
        if (trace) cat("\nForecast GARCH Variance:\n")
        for (i in 1:M) {
            h[N+i] = omega  + sum(beta*h[N+i-(1:q)])
            for (j in 1:p) {
                if (i-j > 0) {
                    s = h[N + i - j]
                } else {
                    s = z[N + i - j]^2
                }
                h[N+i] = h[N+i] + alpha[j] * s
            }
        }
    }
    # Forecast APARCH Variance:
    if (var.model == "aparch") {
        if (trace) cat("\nForecast APARCH Variance:\n")
        for (i in 1:M) {
            h[N+i] = omega  + sum(beta*h[N+i-(1:q)])
            for (j in 1:p) {
                kappa = garchKappa(cond.dist = cond.dist, gamma = gamma[j],
                    delta = delta, skew = skew, shape = shape)
                if (i-j > 0) {
                    s = kappa * h[N + i - j]
                } else {
                    s = (abs(z[N + i - j]) - gamma[j]*z[N + i - j])^delta
                }
                h[N+i] = h[N+i] + alpha[j] * s
            }
        }
    }

    # Forecast and Optionally Trace Mean Model:
    # Note we set maxit=0 to get an object of class Arima with fixed
    #   init parameters ...
    mu <- mu/(1-sum(ar))
    ARMA <- arima(x = object@data, order = c(max(u, 1), 0, max(v, 1)),
                  init = c(ar, ma, mu), transform.pars = FALSE,
                  optim.control = list(maxit = 0))
    prediction = predict(ARMA, n.ahead)
    meanForecast = as.vector(prediction$pred)
    if(mse=="uncond") {
        meanError = as.vector(prediction$se)
    } else {
	# coefficients of h(t+1)
	a_vec <- rep(0,(n.ahead))
	hhat <- h[-(1:N)]^(2/delta[[1]]) #-> [[1]] to omit name of delta
	u2 <- length(ar)
	meanError <- hhat[1]
        a_vec[1] = ar[1] + ma[1]
        meanError <- na.omit(c(meanError,sum(hhat[1:2]*c(a_vec[1]^2,1))))
        if ((n.ahead - 1) > 1) {
            for( i in 2:(n.ahead - 1)) {
                a_vec[i] <- ar[1:min(u2,i-1)]*a_vec[(i-1):(i-u2)] +
                    ifelse(i>u,0,ar[i]) + ifelse(i>v,0,ma[i])
                meanError <- na.omit(c(meanError,
                                       sum(hhat[1:(i+1)]*c(a_vec[i:1]^2,1))))
            }
        }
	meanError <- sqrt(meanError)
    }
    if (trace) {
        cat("\nForecast ARMA Mean:\n")
        print(ARMA)
        cat("\n")
        print(prediction)
    }


    # Standard Deviations:
    standardDeviation = h^(1/delta)

    # Plotting the predictions

    if (plot) {
        if(is.null(nx))
            nx <- round(length(object@data)*.25)
	t <- length(object@data)
	x <- c(object@data[(t-nx+1):t],meanForecast)

	# Computing the appropriate critical values

	if (is.null(conf))
            conf <- 0.95

        if (is.null(crit_val)) {
            if (object@fit$params$cond.dist=="norm") {
                crit_valu <- qnorm(1-(1-conf)/2)
                crit_vald <- qnorm((1-conf)/2)
            }
            if (object@fit$params$cond.dist=="snorm") {
                crit_valu <- qsnorm(1-(1-conf)/2,xi=coef(object)["skew"])
                crit_vald <- qsnorm((1-conf)/2,xi=coef(object)["skew"])
            }
            if (object@fit$params$cond.dist=="ged") {
                crit_valu <- qged(1-(1-conf)/2,nu=coef(object)["shape"])
                crit_vald <- qged((1-conf)/2,nu=coef(object)["shape"])
            }
            if (object@fit$params$cond.dist=="sged") {
                crit_valu <- qsged(1-(1-conf)/2,nu=coef(object)["shape"],
                                   xi=coef(object)["skew"])
                crit_vald <- qsged((1-conf)/2,nu=coef(object)["shape"],
                                   xi=coef(object)["skew"])
            }
            if (object@fit$params$cond.dist=="std") {
                crit_valu <- qstd(1-(1-conf)/2,nu=coef(object)["shape"])
                crit_vald <- qstd((1-conf)/2,nu=coef(object)["shape"])
            }
            if (object@fit$params$cond.dist=="sstd") {
                crit_valu <- qsstd(1-(1-conf)/2,nu=coef(object)["shape"],
                                   xi=coef(object)["skew"])
                crit_vald <- qsstd((1-conf)/2,nu=coef(object)["shape"],
                                   xi=coef(object)["skew"])
            }
            if (object@fit$params$cond.dist=="snig") {
                crit_valu <- qsnig(1-(1-conf)/2,zeta=coef(object)["shape"],
                                   rho=coef(object)["skew"])
                crit_vald <- qsnig((1-conf)/2,zeta=coef(object)["shape"],
                                   rho=coef(object)["skew"])
            }
            if (object@fit$params$cond.dist=="QMLE") {
                e <- sort(object@residuals/object@sigma.t)
                crit_valu <- e[round(t*(1-(1-conf)/2))]
                crit_vald <- e[round(t*(1-conf)/2)]
            }
        } else {
            if (length(crit_val)==2) {
                crit_valu <- crit_val[2]
                crit_vald <- crit_val[1]
            }
            if (length(crit_val)==1) {
                crit_valu <- abs(crit_val)
                crit_vald <- -abs(crit_val)
            }
        }

	int_l <- meanForecast+crit_vald*meanError
	int_u <- meanForecast+crit_valu*meanError
	ylim_l <- min(c(x,int_l)*(.95))
	ylim_u <- max(c(x,int_u)*(1.05))

	plot(x,type='l',ylim=c(ylim_l,ylim_u))
	title("Prediction with confidence intervals")
	lines((nx+1):(nx+n.ahead), meanForecast, col = 2, lwd = 2)
	lines((nx+1):(nx+n.ahead), int_l, col = 3, lwd = 2)
	lines((nx+1):(nx+n.ahead), int_u, col = 4, lwd = 2)
	polygon(c((nx+1):(nx+n.ahead),(nx+n.ahead):(nx+1)),
                c(int_l, int_u[n.ahead:1]),
                border = NA, density = 20, col = 5, angle = 90)
	es1 <- as.expression(substitute(hat(X)[t+h] + crit_valu*sqrt(MSE),
            list(crit_valu=round(crit_valu,3))))
	es2 <- as.expression(substitute(hat(X)[t+h] - crit_vald*sqrt(MSE),
            list(crit_vald=abs(round(crit_vald,3)))))
	es3 <- expression(hat(X)[t+h])
	legend("bottomleft",c(es3,es2,es1),col=2:4,lty=rep(1,3),lwd=rep(2,3))
	grid()
    }


    # Result:
    if(plot) {
	forecast = data.frame(
        meanForecast = meanForecast,
        meanError = meanError,
        standardDeviation = standardDeviation[-(1:N)],
	lowerInterval = int_l,
	upperInterval = int_u)
    } else {
	forecast = data.frame(
        meanForecast = meanForecast,
        meanError = meanError,
        standardDeviation = standardDeviation[-(1:N)])
    }

    # Return Value:
    forecast
})

################################################################################
