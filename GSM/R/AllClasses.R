.onAttach <- function(lib, pkg) {
	packageStartupMessage(sprintf("Package %s (%s) loaded.
To cite, see citation(\"%s\")\n", pkg, packageDescription(pkg)$Version, pkg))
}

setClass(Class = "gsm", representation(fdens = "matrix", theta = "numeric", weight = "matrix", label = "matrix", data = "numeric"))

setMethod("initialize", "gsm",
		function(.Object, fdens = matrix(NA), theta = numeric(0), weight = matrix(NA), label = matrix(NA), data = numeric(0)) {
			.Object@fdens <- fdens
			.Object@theta <- theta
			.Object@weight <- weight
			.Object@label <- label
			.Object@data <- data
			.Object
		}
)

setMethod("summary",
		"gsm",
		function(object, plot = FALSE, start = 1, ...) {
			numsim <- dim(object@fdens)[1]
			if (numsim < start)
				stop(paste("number of draws to use smaller than those available (", numsim, ")", sep = ""))
			y <- object@data
			thetasim <- object@theta[start:numsim]
			weightsim <- object@weight[start:numsim, ]
			if (numsim == start) {
				out <- list(theta = summary(thetasim), `weights posterior means` = weightsim)
			}
			else {
				out <- list(theta = summary(thetasim), `weights posterior means` = colMeans(weightsim))
			}
			if (plot) {
				J <- dim(object@weight)[2]
				if (numsim == start) {
					barplot(weightsim, names.arg = 1:J, xlab = "Mixture component", ylab = "Estimated posterior mean", main = "Mixture weights posterior means")
				}
				else {
					barplot(colMeans(weightsim), names.arg = 1:J, xlab = "Mixture component", ylab = "Estimated posterior mean",
							main = "Mixture weights posterior means")
				}
			}
			out
		}
)

setMethod("predict",
		"gsm",
		function(object, thresh, start = 1, ...) {
			numsim <- dim(object@fdens)[1]
			if (numsim < start)
				stop(paste("number of draws to use smaller than those available (", numsim, ")", sep = ""))
			thetasim <- object@theta[start:numsim]
			weightsim <- object@weight[start:numsim, ]
			runs <- dim(weightsim)[1]
			prob <- vector(length = runs)
			for (i in 1:runs) prob[i] <- 1 - pgsm(thresh, weightsim[i, ], thetasim[i])
			prob
		}
)

setMethod("plot",
		signature(x = "gsm", y = "missing"),
		function(x, ndens = 5, xlab = "x", ylab = "density", nbin = 10, histogram = FALSE, bands = FALSE, confid = .95, start = 1, ...) {
			data <- x@data
			fdens <- x@fdens[-c(1:(start-1)), ]
			numsim <- dim(fdens)[1]
			if (numsim < start)
				stop(paste("number of draws to use smaller than those available (", numsim, ")", sep = ""))
			rpl <- ifelse(ndens > (numsim - start + 1), TRUE, FALSE)
			G <- dim(fdens)[2]
			y.grid <- seq(min(data)*.66, max(data)*1.5, length = G)
			xlim <- range(data)
			ylim <- c(0, max(fdens))
			yval <- apply(fdens, 2, mean)
			if (histogram) {
				hist(data, freq = FALSE, breaks = nbin, col = gray(.5), border = gray(.25), xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, main = "")
			}
			else {
				plot(c(0, y.grid), c(0, yval), type = "n", xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab)
			}
			if (bands) {
				color <- rgb(190, 190, 190, alpha=180, maxColorValue=255)
				grid <- c(y.grid, rev(y.grid))
				perc <- (1 - confid) / 2
				curves <- c(allcurves.q(fdens, perc), rev(allcurves.q(fdens, (1 - perc))))
				polygon(grid, curves, col = color, lty = 1, lwd = 2, border = NA)
			}
			for (i in sample(1:dim(fdens)[1], ndens, replace = rpl)) lines(y.grid, fdens[i, ], col = "red")
			lines(y.grid, yval, lwd = 2)
			rug(data)
			out <- list(xval = y.grid, yval = yval)
			invisible(out)
		}
)
