#' Plot Model Fits
#'
#' Methods to plot growth model fits together with the data and, alternatively,
#' plot diagnostics
#'
#' @param x an object returned by a model fitting function of package
#'   \pkg{growthrates}, that can contain one or multiple fits.
#' @param y (ignored) for compatibility with the default plot method.
#' @param log a character string which contains \code{"y"} if the y axis is to
#'   be logarithmic.
#' @param which either \code{"fit"} (default) or \code{"diagnostics"}.
#' @param \dots other arguments pased to the plotting methods,
#' see \code{\link{plot.default}} and \code{\link{par}}.
#'
#' @details The plot methods detect automatically which type of plot is
#'   appropriate, depending on the class of \code{x} and can plot either one
#'   single model fit or a complete series (multiple fits). In the latter case
#'   it may be wise to redirect the graphics to an external file (e.g. a pdf)
#'   and / or to use tomething like \code{par(mfrow=c(3,3))}.
#'
#'   The \code{lines}-method is currently only available for single fits.
#'
#'   If you need more control, you can of course also write own plotting functions.
#'
#' @examples
#'
#' data(bactgrowth)
#' splitted.data <- multisplit(bactgrowth, c("strain", "conc", "replicate"))
#'
#' ## get table from single experiment
#' dat <- splitted.data[["D:0:1"]]
#'
#' fit1 <- fit_spline(dat$time, dat$value)
#' plot(fit1, log="y")
#' plot(fit1)
#'
#' ## derive start parameters from spline fit
#' p <- coef(fit1)
#'
#' ## subset of first 10 data
#' first10 <-  dat[1:10, ]
#' fit2 <- fit_growthmodel(grow_exponential, p=p, time=first10$time, y=first10$value)
#'
#' p <- c(coef(fit1), K = max(dat$value))
#' fit3 <- fit_growthmodel(grow_logistic, p=p, time=dat$time, y=dat$value, transform="log")
#'
#' plot(fit1)
#' lines(fit2, col="green")
#' lines(fit3, col="red")
#'
#'
#' all.fits <- all_splines(value ~ time | strain + conc + replicate, data = bactgrowth)
#' par(mfrow=c(3,3))
#' plot(all.fits)
#'
#' ## it is also possible to plot a single fit or a subset of the fits
#' par(mfrow=c(1,1))
#' plot(all.fits[["D:0:1"]])
#' par(mfrow=c(2,2))
#' plot(all.fits[1:4])
#'
#' ## plot only the 'R' strain
#' par(mfrow=c(4, 6))
#' plot(all.fits[grep("R:", names(all.fits))])
#'
#' @seealso \code{\link{plot.default}}, \code{\link{par}},
#'   \code{\link{fit_growthmodel}}, \code{\link{fit_easylinear}},
#'   \code{\link{all_growthmodels}}, \code{\link{all_easylinear}}
#' @name plot
#'
NULL

#' @rdname plot
#' @exportMethod plot
#'
setMethod("plot", c("nonlinear_fit", "missing"),
          function(x, y, log="", which=c("fit", "diagnostics"), ...) {
            which <- match.arg(which)

            switch(which,
                   fit = {
                     obs <- obs(x)
                     plot(obs, log=log, ...)
                     times <- seq(min(obs$time), max(obs$time), length.out=200)
                     sim <- x@FUN(times, coef(x))
                     lines(sim[,"time"], sim[,"y"], col="blue", ...)
                   },
                   diagnostics = {
                     opar <- par(no.readonly = TRUE)
                     ## diagnostics from FME: iteration, residuals vs. index
                     plot(x@fit, ..., mfrow=c(2,2))

                     ## residuals vs. fitted
                     obs <- obs(x)
                     sim <- x@FUN(obs$time, coef(x))
                     plot(residuals(x) ~ sim[,"y"], xlab="fitted", ylab="residuals")
                     abline(h=0, col="grey")
                     ## normal q-q-plot
                     qqnorm(residuals(x))
                     qqline(residuals(x))
                     par(opar)
                   }
            )
          }
)


#' @rdname plot
#' @exportMethod lines
#'
setMethod("lines", "nonlinear_fit",
          function(x, ...) {
            obs <- obs(x)
            times <- seq(min(obs$time), max(obs$time), length.out=200)
            sim <- x@FUN(times, coef(x))
            lines(sim[,"time"], sim[,"y"], ...)
          }
)



#' @rdname plot
#'
setMethod("plot", c("easylinear_fit", "missing"),
          function(x, y, log="", which=c("fit", "diagnostics"), ...) {
            which <- match.arg(which)

            switch(which,
                   fit = {
                     obs <- obs(x)
                     plot(obs[,"y"] ~ obs[,"time"], xlab="time", ylab="y",
                          log=log, ...)
                     points(obs[x@ndx,"y"] ~ obs[x@ndx,"time"], pch=16, col="red")
                     time <- seq(min(obs[,"time"]), max(obs[,"time"]), length=200)
                     lines(time, x@FUN(time, coef(x))[,"y"], ...)
                   },
                   diagnostics = {
                     opar <- par(no.readonly = TRUE)
                     on.exit(par(opar))
                     par(mfrow=c(1,2))

                     ## residuals vs. fitted
                     obs <- obs(x)
                     sim <- x@FUN(obs$time, coef(x))
                     plot(residuals(x) ~ fitted(x@fit), xlab="fitted", ylab="residuals")
                     abline(h=0, col="grey")
                     ## normal q-q-plot
                     qqnorm(residuals(x))
                     qqline(residuals(x))
                   }
            )
          }
)


#' @rdname plot
#'
setMethod("plot", c("smooth.spline_fit", "missing"),
        function(x, y, ...) {
            time <- x@obs$time
            y    <- x@obs$y
            px   <- x@xy[1]
            py   <- x@xy[2]
            y0   <- coef(x)["y0"]
            mumax   <- coef(x)["mumax"]

            xnew <- seq(min(time), max(time), length.out = 200)
            ynew <- predict(x@fit, x = xnew)$y
            dspl <- predict(x@fit, x = xnew, deriv=1)

            ## plot results
            plot(time, y, ...) #ylim=c(range(c(y, ynew, dspl$y))), ...)
            lines(xnew, exp(ynew))
            #lines(xnew, dspl$y, col="red")

            ### plot tangent
            points(px, py, pch=16, col="red")
            lines(xnew, y0 * exp(mumax * xnew), col="blue")

            ## draw predictions
            sim <- x@FUN(xnew, coef(x))
            lines(sim[,"time"], sim[,"y"], col="blue")

            ## alternative for log-transformed scale
            #xx <- px + c(-4, 4)
            #yy <- log(py) + x@par["mumax"] * c(-4, 4)
            #lines(xx, exp(yy), col="grey")

            invisible()
          }


)


#' @rdname plot
#'
setMethod("lines", "easylinear_fit",
          function(x, ...) {
            obs <- obs(x)
            ## mark used data points
            points(obs[x@ndx,"y"] ~ obs[x@ndx,"time"], ...)

            ## draw predictions
            time <- seq(min(obs[,"time"]), max(obs[,"time"]), length=200)
            lines(time, x@FUN(time, coef(x))[,"y"], ...)
          }
)



### ----------------------------------------------------------------------------
### Plotting of Multiple Fits
### ----------------------------------------------------------------------------


#' @rdname plot
#'
setMethod("plot", c("multiple_fits", "missing"),
          function(x, y, ...) {
            ## transfer both, list element and name to plot
            lapply(seq_along(x@fits),
                   FUN = function(i, y, n)
                     plot(y@fits[[i]], main = n[i], ...),
                   y=x, n=names(x@fits))
            invisible()
          })



#setMethod("plot", "multiple_easylinear_fits",
#          function(x, ...) {
#            lapply(x@fits, plot, ...)
#            invisible()
#          })




