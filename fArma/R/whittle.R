
## # This program is free software; you can redistribute it and/or modify
## # it under the terms of the GNU General Public License as published by
## # the Free Software Foundation; either version 2, or (at your option)
## # any later version.
## #
## # This program is distributed in the hope that it will be useful, but
## # WITHOUT ANY WARRANTY; without even the implied warranty of
## # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## # General Public License for more details.
## #
## # A copy of the GNU General Public License is available via WWW at
## # http://www.gnu.org/copyleft/gpl.html.  You can also obtain it by
## # writing to the Free Software Foundation, Inc., 59 Temple Place,
## # Suite 330, Boston, MA  02111-1307  USA.

## # Copyrights (C)
## # for this R-port:
## #   1999 - 2007, Diethelm Wuertz, GPL
## #   Diethelm Wuertz <wuertz@itp.phys.ethz.ch>
## #   info@rmetrics.org
## #   www.rmetrics.org
## # for the code accessed (or partly included) from other R-ports:
## #   see R's copyright and license files
## # for the code accessed (or partly included) from contributed R-ports
## # and other sources
## #   see Rmetrics's copyright file


## ################################################################################
## # FUNCTIONS:            WHITTLE ESTIMATOR:
## #  whittleFit            Whittle Estimator
## #  .CetaFGN               Internal Functions ...
## #  .CetaARIMA
## #  .Qeta
## #  .fspecFGN
## #  .ffourier.FGN.est
## #  .FGN.spectrum
## #  .FGN.B.est.adjust
## #  .FGN.B.est
## #  .fspecARIMA
## #  .per
## #  .Qmin2
## #  .whittle
## ################################################################################


## ################################################################################
## # DESCRIPTION:
## #   The functions are reimplemented from the appendix of J. Beran "Statistics
## #   for long-memory processes", Chapman and Hall 1984
## # LICENSE:
## #   Permission is hereby given to StatLib to redistribute this software.
## #   The software can be freely used for non-commercial purposes, and can
## #   be freely distributed for non-commercial purposes only.
## # AUTHORS:
## #   Jan Beran <jberan@iris.rz.uni-konstanz.de>
## #   Modified: Martin Maechler <maechler@stat.math.ethz.ch>
## #   Modified: Diethelm Wuertz <wuertz@itp.phys.ethz.ch> for this R-Port


## whittleFit =
## function(x, order = c(1, 1), subseries = 1, method = c("fgn", "farma"),
## trace = FALSE, spec = FALSE, title = NULL, description = NULL)
## {   # A function implemented by Diethelm Wuertz

##     # Description:
##     #   Minimizes an approximate log-likelihood function applied to the
##     #   spectral density to obtain an estimate of the parameters of a
##     #   process.

##     # Details:
##     #   Function and programs for the calculation of Whittle's estimator
##     #   and the goodness of fit statistic as defined in Beran (1992). The
##     #   models are fractional Gaussian noise or fractional Arima. The data
##     #   series may be divided into subseries for which the parameters are
##     #   fitted separately.
##     #
##     #   There are several options for using the Whittle estimator. Some are
##     #   described below.
##     #   1.  One can optionally subdivide the series into "subseries".
##     #   3.  One can output the periodogram by "spec".
##     #   4.  One can "trace" intermediate minimization results.
##     #   5.  The "model" can be either farma or fgn.
##     #   6.  If the model is farma, the "order" has to be specified.
##     #   7.  The starting value of H for the minimization procedure is "h".
##     #   8.  "ar" and "ma" are starting values of the time series coefficients.
##     #       (Length of vectors should be the same as p and q).

##     # FUNCTION:

##     # Settings:
##     data = list(x = x)
##     x = as.vector(x)

##     # Start Values:
##     h = 0.7
##     ar = rep(0.5, length = order[1]) / order[1]
##     ma = rep(0.5, length = order[2]) / order[2]

##     # Estimate:
##     if(trace) cat("Iteration Path:\n")
##     result = .whittle(xinput = x, nsub = subseries, model = method[1],
##         pp = order[1], qq = order[2], h = h, ar = ar, ma = ma, out = trace,
##         spec = spec)[[1]]
##     result$H = result$par
##     result$par = NULL

##     # Add:
##     if(is.null(title)) title = "Hurst Exponent from Whittle Estimator"
##     if(is.null(description)) description = description()

##     # Return Value:
##     new("fHURST",
##         call = match.call(),
##         method = paste(method[1], "whittle"),
##         hurst = result,
##         parameter = list(subseries = subseries, order = order,
##             h = h, ar = ar, ma = ma),
##         data = data,
##         fit = result,
##         plot = list(doplot = FALSE),
##         title = title,
##         description = description
##         )
## }


## ################################################################################
## # Functions to make this function independent from Beran's code


## .CetaFGN =
## function(eta)
## {   # A function implemented by Diethelm Wuertz

##     # Description:
##     #   Computes covariance matrix of hat{eta} for fGn

##     # Author:
##     #   Jan Beran; modified: Martin Maechler Sep 95, Diethelm Wuertz 2004

##     # FUNCTION:

##     # Settings:
##     M = length(eta)

##     # Size of steps in Riemann sum: 2*pi/m
##     m = 10000
##     mhalfm = trunc((m-1)/2)

##     # Size of delta for numerical calculation of derivative
##     delta = 1.0e-9

##     # Partial derivatives of log f (at each Fourier frequency)
##     lf = matrix(1, ncol = M, nrow = mhalfm)
##     f0 = .fspecFGN(eta,m)$fspec
##     for (j in (1:M)) {
##         etaj = eta
##         etaj[j] = etaj[j] + delta
##         fj = .fspecFGN(etaj, m)$fspec
##         lf[,j] = log(fj/f0)/delta }

##         # Calculate D:
##     Djl = matrix(1,ncol = M, nrow = M)
##     for (j in (1:M)) {
##         for(l in (1:M)) {
##             Djl[j,l] = 2*2*pi/m*sum(lf[,j]*lf[,l])
##         }
##     }
##     ans = drop(matrix(4*pi*solve(Djl), ncol = M, nrow = M, byrow = TRUE))

##     # Return Value:
##     ans
## }


## # ------------------------------------------------------------------------------


## .CetaARIMA =
## function(eta, p, q)
## {   # A function implemented by Diethelm Wuertz

##     # Description:
##     #   Computes ovariance matrix of hat{eta} for fractional ARIMA

##     # Author:
##     #   Jan Beran; modified: Martin Maechler, Diethelm Wuertz 2004

##     # FUNCTION:

##     # Settings:
##     M = length(eta)

##     # Size of steps in Riemann sum: 2*pi/m:
##     m = 10000
##     mhalfm = trunc((m-1)/2)

##     # Size of delta for numerical calculation of derivative:
##     delta = 1.0e-9
##     # partial derivatives of log f (at each Fourier frequency)
##     lf = matrix(1, ncol = M, nrow = mhalfm)
##     f0 = .fspecARIMA(eta, p, q, m)$fspec
##     for (j in (1:M)) {
##         etaj = eta
##         etaj[j] = etaj[j]+delta
##         fj = .fspecARIMA(etaj, p, q, m)$fspec
##         lf[,j] = log(fj/f0)/delta
##     }

##     # Calculate D:
##     Djl = matrix(1,ncol = M, nrow = M)
##     for (j in (1:M)) {
##         for (l in (1:M)) {
##             Djl[j,l] = 2*2*pi/m*sum(lf[,j]*lf[,l])
##         }
##     }
##     ans = drop(matrix(4*pi*solve(Djl),ncol = M, nrow = M, byrow = TRUE))

##     # Return Value:
##     ans
## }


## # ------------------------------------------------------------------------------


## .Qeta =
## function(eta)
## {   # A function implemented by Diethelm Wuertz

##     # Description:
##     #   Calculation of A, B and Tn = A/B**2
##     #   where A = 2pi/n sum 2*[I(lambda = j)/f(lambda = j)],
##     #         B = 2pi/n sum 2*[I(lambda = j)/f(lambda = j)]**2  and
##     #   the sum is taken over all Fourier frequencies
##     #   lambda = j = 2pi*j/n (j=1,...,(n-1)/2.
##     #   f is the spectral density of fractional Gaussian
##     #   noise or fractional ARIMA(p,d,q) with self-similarity parameter H = h.
##     #   cov(X(t),X(t+k))=integral(exp(iuk)f(u)du)

##     # Arguments:
##     #   h
##     #   (n, nhalfm = trunc[(n-1)/2] and the
##     #   nhalfm-dimensional  GLOBAL vector `yper' must be defined.)

##     # Value:
##     #   list(n=n,h=h,A=A,B=B,Tn=Tn,z=z,pval=pval, theta1=theta1,fspec=fspec)
##     #   Tn is the goodness of fit test statistic
##     #   Tn=A/B**2 defined in Beran (1992),
##     #   z is the standardized test statistic,
##     #   pval the corresponding p-value P(w>z).
##     #   theta1 is the scale parameter such that
##     #   f=theta1*fspec and integral(log[fspec]) = 0.

##     # Note:
##     #   yper[1] must be the periodogram I(lambda = 1) at
##     #   the frequency 2pi/n (i.e. not the frequency zero !).

##     # Author:
##     #   Jan Beran; modified: Martin Maechler Sep. 95, Diethelm Wuertz 2004

##     # FUNCTION:

##     # To Suppress No Visible Bindings Warning/Error:
##     if(FALSE) { imodel = n = p = yper = NA }

##     # Settings:
##     h = eta[1]
##     if(imodel == 1) {
##         fspec = .fspecFGN(eta, n)
##         theta1 = fspec$theta1
##         fspec = fspec$fspec
##     } else {
##         fspec = .fspecARIMA(eta, p, q, n)
##         theta1 = fspec$theta1
##         fspec = fspec$fspec
##     }
##     yf = yper/fspec
##     yfyf = yf**2
##     A = 2*(2*pi/n)*sum(yfyf)
##     B = 2*(2*pi/n)*sum(yf)
##     Tn = A/(B**2)
##     z = sqrt(n)*(pi*Tn-1)/sqrt(2)
##     pval = 1-pnorm(z)
##     theta1 = B/(2*pi)
##     fspec = fspec
##     Qresult = list(n = n, h = h, eta = eta, A = A,B = B, Tn = Tn,
##         z = z, pval = pval, theta1 = theta1, fspec = fspec)
##     ans = drop(Qresult)

##     # Return value:
##     ans
## }


## # ------------------------------------------------------------------------------


## .fspecFGN =
## function(eta, m)
## {   # A function implemented by Diethelm Wuertz

##     # Description:
##     #   Calculation of the spectral density f of normalized fractional
##     #   Gaussian noise with self-similarity parameter H=h at the
##     #   Fourier frequencies 2*pi*j/m (j=1,...,(m-1)).

##     # Arguments:
##     #   m = sample size
##     #   h = self-similarity parameter

##     # Value:
##     #   list(fspec = fspec, theta1 = theta1)

##     # Note:
##     #   1. cov(X(t),X(t+k)) = integral[exp(iuk)f(u)du]
##     #   2. f = theta1*fspec and integral[log(fspec)] = 0.

##     # Author:
##     #   Jan Beran; modified: Martin Maechler, Diethelm Wuertz 2004

##     # FUNCTION:

##     # Settings:

##     # Taqqu: This Implementation is more efficient than Beran's:
##     fspec = .ffourier.FGN.est(eta, m)
##     logfspec = log(fspec)
##     fint = 2/(m)*sum(logfspec)
##     theta1 = exp(fint)
##     fspec = fspec/theta1
##     ans = drop(list(fspec = fspec, theta1 = theta1))

##     # Return Value:
##     ans
## }


## # ------------------------------------------------------------------------------


## .ffourier.FGN.est =
## function(H, n)
## {   # A function implemented by Diethelm Wuertz

##     # Description:
##     #   Internal Function

##     # Author:
##     #   Jan Beran; modified: Martin Maechler, Diethelm Wuertz 2004

##     # FUNCTION:

##     # Spectrum:
##     ans = .FGN.spectrum((2 * pi * (1:((n - 1)/2)))/n, H)/pi/2

##     # Return Value:
##     ans
## }


## # ------------------------------------------------------------------------------


## .FGN.spectrum =
## function(lambda, H)
## {   # A function implemented by Diethelm Wuertz

##     # Description:
##     #   Internal Function

##     # Author:
##     #   Jan Beran; modified: Martin Maechler, Diethelm Wuertz 2004

##     # FUNCTION:

##     # Settings:
##     ans = 2*sin(pi*H)*gamma(2*H+1)*(1-cos(lambda))*(lambda^(-2*H-1) +
##         .FGN.B.est.adjust(lambda, H))
##     ans
## }


## # ------------------------------------------------------------------------------


## .FGN.B.est.adjust =
## function(lambda, H)
## {   # A function implemented by Diethelm Wuertz

##     # Description:
##     #   Internal Function

##     # Author:
##     #   Jan Beran; modified: Martin Maechler, Diethelm Wuertz 2004

##     # FUNCTION:

##     # Settings:
##     B = .FGN.B.est(lambda, H)
##     ans = (1.0002-0.000134*lambda) * (B-2^(-7.65*H-7.4))
##     ans
## }


## # ------------------------------------------------------------------------------


## .FGN.B.est =
## function(lambda, H)
## {   # A function implemented by Diethelm Wuertz

##     # Description:
##     #   Internal Function

##     # Author:
##     #   Jan Beran; modified: Martin Maechler, Diethelm Wuertz 2004

##     # FUNCTION:

##     # Settings:
##     d = -(2*H+1)
##     dprime = -2*H
##     a = function(lambda, k) { 2 * k * pi+lambda }
##         b = function(lambda, k) { 2 * k * pi-lambda }
##         a1 = a(lambda, 1); b1 = b(lambda, 1)
##         a2 = a(lambda, 2); b2 = b(lambda, 2)
##     a3 = a(lambda, 3)
##     b3 = b(lambda, 3)
##     a4 = a(lambda, 4)
##     b4 = b(lambda, 4)
##     ans = a1^d+b1^d+a2^d+b2^d+a3^d+b3^d+
##         (a3^dprime+b3^dprime+a4^dprime+b4^dprime)/(8*pi*H)
##     ans
## }


## # ------------------------------------------------------------------------------


## .fspecARIMA =
## function(eta, p, q, m)
## {   # A function implemented by Diethelm Wuertz

##     # Description:
##     #   Internal Function

##     # Author:
##     #   Jan Beran; modified: Martin Maechler, Diethelm Wuertz 2004

##     # FUNCTION:

##     # Settings:
##     h = eta[1]
##     phi = c()
##     psi = c()
##     mhalfm = trunc((m-1)/2)
##     x = 2*pi/m*(1:mhalfm)
##     # Calculation of f at Fourier frequencies
##     far = (1:mhalfm)/(1:mhalfm)
##     fma = (1:mhalfm)/(1:mhalfm)
##     if(p > 0) {
##         phi = cbind(eta[2:(p+1)])
##         cosar = cos(cbind(x) %*% rbind(1:p))
##         sinar = sin(cbind(x) %*% rbind(1:p))
##         Rar = cosar %*% phi
##         Iar = sinar %*% phi
##         far = (1-Rar)**2 + Iar**2 }
##     if(q > 0) {
##         psi = cbind(eta[(p+2):(p+q+1)])
##         cosar = cos(cbind(x) %*% rbind(1:q))
##         sinar = sin(cbind(x) %*% rbind(1:q))
##         Rar = cosar %*% psi
##         Iar = sinar %*% psi
##         fma = (1+Rar)**2 + Iar**2 }
##     fspec = fma/far*sqrt((1-cos(x))**2 + sin(x)**2)**(1-2*h)
##     theta1 = 1/(2*pi)
##     ans = list(fspec = fspec, theta1 = theta1)
##     ans
## }


## # ------------------------------------------------------------------------------


## .per =
## function(z)
## {   # A function implemented by Diethelm Wuertz

##     # Description:
##     #   Internal Function

##     # Author:
##     #   Jan Beran; modified: Martin Maechler, Diethelm Wuertz 2004

##     # FUNCTION:

##     # Settings:
##     n = length(z)
##     ans = (Mod(fft(z))**2/(2*pi*n))[1:(n %/% 2 + 1)]

##     # Return Value:
##     ans
## }


## # ------------------------------------------------------------------------------


## .Qmin2 =
## function(etatry)
## {   # A function implemented by Diethelm Wuertz

##     # Description:
##     #   Internal Function

##     # FUNCTION:

##     # Compute:
##     ans = .Qeta(etatry)$B
##     assign("bBb", ans, pos = 1)

##     # Return Value:
##     ans
## }


## # ------------------------------------------------------------------------------
## # Internal Function: WHITTLE ESTIMATOR:


## .whittle =
## function(xinput, nsub = 1, model = c("farma", "fgn"), pp = 1,
## qq = 1, h = 0.5, ar = c(0.5), ma = c(0.5), out = TRUE, spec = FALSE)
## {   # A function implemented by Diethelm Wuertz

##     # Description:
##     #   Internal Function

##     # Author:
##     #   Jan Beran; modified: Martin Maechler, Diethelm Wuertz 2004

##     # FUNCTION:

##     # To suppress No Visible Bindings Warning/Error:
##     if(FALSE) { n = p = q = imodel = out = yper = bBb = NA }

##     # Settings:
##     model = model[1]
##     assign("out", out, pos = 1)
##     nmax = length(xinput)
##     startend = c(1, nmax)
##     istart = startend[1]
##     iend = startend[2]
##     nloop = nsub
##     assign("n", trunc((iend - istart+1)/nloop), pos = 1)
##     nhalfm = trunc((n - 1)/2)
##     if(model == "farma")
##         assign("imodel", 2, pos = 1)
##     if(model == "fgn")
##         assign("imodel", 1, pos = 1)
##     assign("p", 0, pos = 1)
##     assign("q", 0, pos = 1)
##     p <- q <- 0
##     if(imodel == 2) {
##         assign("p", pp, pos = 1); assign("q", qq, pos = 1)
##     } else {
##         assign("p", 0, pos=1); assign("q", 0, pos = 1)
##     }
##     eta = c(h)
##     if(p > 0) eta[2:(p+1)] = ar
##     if(q > 0) eta[(p+2):(p+q+1)] = ma
##     M = length(eta) #loop
##     thetavector = c()
##     i0 = istart
##     flax = vector("list", nloop)

##     # Necessary to make nsub/nloop work.  VT.
##     for (iloop in (1:nloop)) {

##         h = max(0.2, min(h, 0.9))
##         eta[1] = h
##         i1 = i0+n - 1
##         y = xinput[i0:i1]

##         # Standardize Data:
##         vary = var(y)
##         y = (y - mean(y))/sqrt(var(y))

##         # Periodogram of the Data:
##         if(spec) {
##             assign("yper", .per(y)[2:(nhalfm+1)], pos = 1)
##         } else {
##             assign("yper", .per(y)[2:(nhalfm+1)], pos = 1)
##         }
##         s = 2*(1-h)
##         etatry = eta

##         # Modified to make optim not give incorrect result.  VT
##         if(imodel == 1) {
##             result = optim(par = etatry, fn = .Qmin2,
##                 method = "L-BFGS-B", lower = 0, upper = 0.999)
##         } else {
##             result = optim(par = etatry, fn = .Qmin2)
##         }
##         eta = result$par
##         sturno = result$message
##         theta1 = .Qeta(eta)$theta1
##         theta = c(theta1, eta)
##         thetavector = c(thetavector, theta)

##         # Calculate goodness of fit statistic
##         Qresult = .Qeta(eta)    #output
##         M = length(eta)
##         if(imodel == 1) {
##             SD = .CetaFGN(eta)
##             SD = matrix(SD, ncol = M, nrow = M, byrow = TRUE) / n
##         } else {
##             # Changed to eliminate crashing in solve.qr in CetaARIMA.  VT
##             cat("M =", M, "\n")
##             if(M > 2) {
##                 for (i in 3:M) {
##                     for (j in 2:(i-1)) {
##                         temp = eta[i]+eta[j]
##                         if(abs(temp) < 0.0001) {
##                             cat("Problem with estimating confidence intervals,",
##                                 "parameter ", i, "and  parameter ", j,
##                                 "are the same, eliminating.\n")
##                             eta = eta[ - i]
##                             eta = eta[ - j]
##                             M = M - 2
##                             p = p - 1
##                             q = q - 1
##                         }
##                     }
##                 }
##             }
##             SD = .CetaARIMA(eta, p, q)
##             SD = matrix(SD, ncol = M, nrow = M, byrow = TRUE)/n
##         }
##         Hlow = eta[1] - 1.96 * sqrt(SD[1, 1])
##         Hup = eta[1]+1.96 * sqrt(SD[1, 1])
##         if(out) {
##             cat("theta =", theta, fill = TRUE)
##             cat("H =", eta[1], fill = TRUE)
##             cat("95%-CI for H: [", Hlow, ",", Hup, "]", fill = TRUE)
##         }

##         # Changing of the signs of the moving average parameters
##         # in order to respect the sign of the Splus convention
##         if(q > 0) eta[(p+2):(p+q+1)] = -eta[(p+2):(p+q+1)]
##         etalow = c()
##         etaup = c()
##         for (i in (1:length(eta))) {
##             etalow = c(etalow, eta[i] - 1.96 * sqrt(SD[i, i]))
##             etaup = c(etaup, eta[i]+1.96 * sqrt(SD[i, i]))
##         }
##         if(out) {
##             cat("95%-CI:", fill = TRUE)
##             print(cbind(etalow, etaup), fill = TRUE)
##         }
##         if(spec) {
##             cat("Periodogram is in yper", fill = TRUE)
##             assign("fest", Qresult$theta1 * Qresult$fspec, pos = 1)
##             cat("Spectral density is in fest", fill = TRUE)
##         }
##         flax[[iloop]] = list()
##         flax[[iloop]]$par = eta
##         flax[[iloop]]$sigma2 = bBb * var(xinput)
##         flax[[iloop]]$conv.type = sturno
##         remove("bBb", pos = 1)

##         # Next subseries:
##         i0 = i0+n

##         # Changing of the signs of the moving average parameters:
##         if(q > 0) eta[(p+2):(p+q+1)] = -eta[(p+2):(p+q+1)]
##     } # end of nloop


##     # Return:
##     return(flax)

##     # Clean up:
##     remove("n", pos = 1)
##     remove("p", pos = 1)
##     remove("q", pos = 1)
##     remove("imodel", pos = 1)
##     remove("out", pos = 1)
##     if(spec == FALSE) remove("yper", pos = 1)
## }


## ################################################################################
