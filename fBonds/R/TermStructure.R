
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
#   1999 - 2006, Diethelm Wuertz, GPL
#   Diethelm Wuertz <wuertz@itp.phys.ethz.ch>
#   info@rmetrics.org
#   www.rmetrics.org
# for the code accessed (or partly included) from other R-ports:
#   see R's copyright and license files
# for the code accessed (or partly included) from contributed R-ports
# and other sources
#   see Rmetrics's copyright file


################################################################################
# FUNCTION:                     DESCRIPTION:
#   NelsonSiegel                 Nelson Siegel Term Structure
#   Svensson                     Nelson-Siegel-Svensson Term Structure
################################################################################


NelsonSiegel =
function(rate, maturity, doplot = TRUE)
{   # A function written by Diethelm Wuertz

    # Description:
    #   Fit the Yield Curve by the Nelson-Siegel Method

    # Details:
    #   This function finds a global solution. The betas are solved
    #   exactly as a function of tau.

    # Notes:
    #   Prelimiary Status

    # FUNCTION:

    # Settings:
    # Start function minimum (fmin), global minimum (gmin), tau
    # values, and delta tau (dtau)
    n = length(maturity)
    fmin = matrix(rep(NA, times = n), byrow = TRUE, nrow = n)
    gmin = 1.0e99
    tau = rep(NA, times = n)

    # Run the loop over the grid - This gives the start solution
    for (i in 1:n) {
        tau[i] = maturity[i]
        x = maturity/tau[i]
        a = matrix(rep(NA, times = 9), byrow = TRUE, nrow = 3)
            a[1,1] = 1
            a[1,2] = a[2,1] = mean( (1-exp(-x))/x )
            a[1,3] = a[3,1] = mean( (1-exp(-x))/x - exp(-x) )
            a[2,2] = mean( ((1-exp(-x))/x)^2 )
            a[2,3] = a[3,2] = mean( ((1-exp(-x))/x ) * ((1-exp(-x))/x -
                exp(-x)) )
            a[3,3] = mean( ((1-exp(-x))/x - exp(-x))^2 )
        b = c(
            mean ( rate ),
            mean ( rate *  ((1-exp(-x))/x) ),
            mean ( rate * (((1-exp(-x))/x - exp(-x))) ) )
        beta = solve(a, b)
        yfit = beta[1] + beta[2]*exp(-x) + beta[3]*x*exp(-x)
        fmin[i] = sum( (rate - yfit)^2 )
        if (fmin[i] < gmin) {
            gmin = fmin[i]
            gvec = c(beta, tau[i])
        }
    }

    # If desired plot OLS(tau1, tau2):
    if (doplot) {
        plot(tau, log(fmin), type = "b", main = "OLS: Nelson-Siegel")
    }

    # Internal Functions:
    fx = function(maturity, x) {
        x[1] +
        x[2] *  (1-exp(-maturity/x[4]))/(maturity/x[4]) +
        x[3] * ((1-exp(-maturity/x[4]))/(maturity/x[4]) - exp(-maturity/x[4]))
    }
    func = function(x, rate, maturity) {
        sum( (rate -
            x[1] -
            x[2] *  (1-exp(-maturity/x[4]))/(maturity/x[4]) -
            x[3] * ((1-exp(-maturity/x[4]))/(maturity/x[4]) -
                exp(-maturity/x[4]) ) )^2)
    }

    # Optimization:
    # DW: gvec = c(0.0840, -0.0063, 0.0044, 1.7603)
    fit = nlminb(start = gvec, objective = func,
        rate = rate, maturity = maturity,
        control = list(eval.max = 2000, iter.max = 2000))
    names(fit$par) = c("beta0", "beta1", "beta2", "tau1")

    # Plot Curve if desired:
    if (doplot) {
        yfit = fx(maturity, fit$par)
        plot(maturity, rate,
            ylim = c(min(c(rate, yfit)), max(c(rate, yfit)) ),
            main = "Nelson-Siegel" )
        lines(maturity, yfit, col = "steelblue")
        grid()
    }

    # Return Value:
    fit
}


# ------------------------------------------------------------------------------


Svensson =
function(rate, maturity, doplot = TRUE)
{   # A function written by Diethelm Wuertz

    # Description:
    #   Fits the Yield Curve by the Nelson-Siegel-Svensson Method
    # Details:
    #   This function finds a global solution. The betas are solved
    #   exactly as a function of the taus.

    # Notes:
    #   Prelimiary Status

    # FUNCTION:

    # Settings
    # Start function minimum (fmin), global minimum (gmin), tau
    # values, and delta tau (dtau)
    n = length(maturity)
    gmin = 1e99
    fmin = matrix(rep(gmin, times = n*n), byrow = TRUE, nrow = n)
    beta0 = beta1 = beta2 = beta3 = fmin2 = fmin
    tau1 = tau2 = rep(0, times = n)

    # Run the loops over the grid - This gives the start solution
    for (i in 1:n) {
        tau1[i] = maturity[i]
        x1 = maturity/tau1[i]
        for (j in 1:n) {
            tau2[j] = maturity[j]
            x2 = maturity/tau2[j]
            if (i == j) {
                a = matrix(rep(NA, times = 9), byrow = TRUE, nrow = 3)
                    a[1,1] = 1
                    a[1,2] = a[2,1] = mean(exp(-x1))
                    a[1,3] = a[3,1] = mean(2*x1*exp(-x1))
                    a[2,2] = mean(exp(-2*x1))
                    a[2,3] = a[3,2] = mean(2*x1*exp(-2*x1))
                    a[3,3] = mean(4*x1*x1*exp(-2*x1))
                b = c(mean(rate), mean(rate*exp(-x1)),
                    mean(2*rate*x1*exp(-x1)))
                beta = solve(a, b)
                yfit = beta[1] + beta[2]*exp(-x1) + 2*beta[3]*x1*exp(-x1)
                # fmin[i,j] = sum((rate-yfit)^2)
                fmin[i,j] = sum(abs((rate-yfit)))
                beta0[i,j] = beta[1]
                beta1[i,j] = beta[2]
                ## DW: beta2[i,j] = beta[3]
                ## DW: beta3[i,j] = beta[3]
                beta2[i,j] = beta[3]/2
                beta3[i,j] = beta[3]/2
                if (fmin[i,j] < gmin) {
                    gmin = fmin[i,j]
                    gvec = c(beta, beta[length(beta)], tau1[i], tau2[j])
                    gindex = c(i, j)
                }
            } else {
                a = matrix(rep(NA, times = 16), byrow = TRUE, nrow = 4)
                    a[1,1] = 1
                    a[1,2] = a[2,1] = mean(exp(-x1))
                    a[1,3] = a[3,1] = mean(x1*exp(-x1))
                    a[1,4] = a[4,1] = mean(x2*exp(-x2))
                    a[2,2] = mean(exp(-2*x1))
                    a[2,3] = a[3,2] = mean(x1*exp(-2*x1))
                    a[2,4] = a[4,2] = mean(x2*exp(-x1-x2))
                    a[3,3] = mean(x1*x1*exp(-2*x1))
                    a[3,4] = a[4,3] = mean(x1*x2*exp(-x1-x2))
                    a[4,4] = mean(x2*x2*exp(-2*x2))
                b = c(mean(rate), mean(rate*exp(-x1)),
                    mean(rate*x1*exp(-x1)), mean(rate*x2*exp(-x2)))
                beta = solve(a, b)
                yfit = beta[1] + beta[2]*exp(-x1) + beta[3]*x1*exp(-x1) +
                    beta[4]*x2*exp(-x2)
                # fmin[i,j] = sum((rate-yfit)^2)
                fmin[i,j] = sum(abs((rate-yfit)))
                beta0[i,j] = beta[1]
                beta1[i,j] = beta[2]
                beta2[i,j] = beta[3]
                beta3[i,j] = beta[4]
                if (fmin[i,j] < gmin) {
                    gmin = fmin[i,j]
                    gvec = c(beta, tau1[i], tau2[j])
                    gindex = c(i, j)
                }
            }
        }
    }

    # If desired plot "smoothed" OLS(tau1, tau2):
    ntau = length(tau1)-1
    fmin2 = fmin
    # Smoothing:
    for (i in 2:ntau)
        for (j in 2:ntau)
            fmin2[i,j] = ( 4 * fmin[i,j] +
                2 * (fmin[i+1,j] + fmin[i-1,j] + fmin[i,j+1] + fmin[i,j-1] ) +
                (fmin[i+1,j+1] + fmin[i-1,j-1] + fmin[i-1,j+1] + fmin[i+1,j-1] )
                ) / 16
    fmin = fmin2
    if (doplot) {
        tau = 1:length(tau1)
        # tau = tau1
        Z = log(fmin)
        # for (i in tau) Z[i, i] = NA
        # Perspective Plot:
        persp(tau, tau, Z,
            theta = -70, phi = 50,
            xlab = "tau 1", ylab = "tau 2", zlab = "log(fmin)",
            col = "steelblue")
        title(main = "OLS: Nelson-Siegel-Svensson")
        # Image/Contour Plot:
        image(tau, tau, Z, # xlim = c(15, 30), ylim = c(30, 48),
            xlab = "tau 1", ylab = "tau 2")
        contour(tau, tau, log(fmin), nlevels = 50,
            xlab = "tau 1", ylab = "tau 2", add = TRUE)
        points(tau[gindex[1]], tau[gindex[2]], pch = 19, cex = 2,
            col = "steelblue")
        title(main = "OLS: Nelson-Siegel-Svensson")
    }

    # Function to optimize:
    fx = function(maturity, x) {
        x[1] + x[2]*exp(-maturity/x[5]) +
        x[3]*(maturity/x[5])*exp(-maturity/x[5]) +
        x[4]*(maturity/x[6])*exp(-maturity/x[6])
    }
    func = function(x, rate, maturity) {
        sum((rate - x[1] - x[2]*exp(-maturity/x[5]) -
            x[3]*(maturity/x[5])*exp(-maturity/x[5]) -
            x[4]*(maturity/x[6])*exp(-maturity/x[6]))^2)
        sum(abs((rate - x[1] - x[2]*exp(-maturity/x[5]) -
            x[3]*(maturity/x[5])*exp(-maturity/x[5]) -
            x[4]*(maturity/x[6])*exp(-maturity/x[6]))))
    }

    # Optimization:
    fit = nlminb(start = gvec, objective = func,
        rate = rate, maturity = maturity,
        control = list(eval.max = 2000, iter.max = 2000))
    names(fit$par) = c("beta0", "beta1", "beta2", "beta3", "tau1", "tau2")

    # Plot Curve if desired:
    if (doplot) {
        yfit = fx(maturity, fit$par)
        plot(maturity, rate,
            ylim = c(min(c(rate, yfit)), max(c(rate, yfit))),
            main = "Nelson-Siegel-Svensson" )
        grid()
        lines(maturity, yfit, col = "steelblue")
    }

    # Return Value:
    fit
}


################################################################################

