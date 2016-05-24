occ.ref <- function (VT, n, k)
{
	# Fit #
	VND <- matrix(VT[1, ], ncol = 1 + k)
	VT_minus_VND <- VT - matrix(rep(VND, n), byrow = TRUE, ncol = 1 + k)
	VS <- matrix(rep(VT_minus_VND[, 1], 1 + k), ncol = 1 + k)
	coefficients <- 1 - VT_minus_VND / VS
	coefficients[1, ] <- 0

	# Format results #
	x <- list(
		VT = VT,
		coefficients = coefficients, 
		VND = VND,
		VS = VS,
		sigma = matrix(rep(0, 1 + k), ncol = 1 + k)
	)
	x
}

occ.ols <- function (VT, n, k)
{
	# Fit #
	coefficients <- VS <- matrix(nrow = n, ncol = 1 + k)
	coefficients[, 1] <- 0
	VND <- sigma <- matrix(ncol = 1 + k)
	sigma[1, 1] <- 0
	for (i in 1:k) {
		m <- lm(VT[, i + 1] ~ VT[, 1])
		mOcc <- 1 - m$coefficients[2]
		mVND <- m$coefficients[1] / mOcc
		coefficients[, i + 1] <- mOcc
		VND[1, i + 1] <- mVND
		VS[, i + 1] <- (m$fitted.values - mVND) / (1 - mOcc)
		sigma[1, i + 1] <- sqrt(sum(m$residuals^2) / (n - 2))
	}

	# Format results #
	x <- list(
		VT = VT,
		coefficients = coefficients, 
		VND = VND,
		VS = VS,
		sigma = sigma
	)
	if (k == 1) {
		x$VND[1] <- x$VND[2]
		x$VS[, 1] <- x$VS[, 2]
	}
	x
}

occ.reml <- function (VT, n, k)
{
	# Fit #
	f <- nlm(
		function (x)
		{
			# Get parameters #
			Occ <- x[1:k]
			if (min(Occ) < 0 | max(Occ) > 1) {
				return(9e+30)
			}
			VND <- x[k + 1]
			sigma2 <- x[k + 2]^2
			VS <- x[(k + 3):(k + 3 + n - 1)]
			if (min(VS) < 0) {
				return(9e+30)
			}
			CO <- c(1, 1 - Occ)
			mVT <- apply(VT, 2, mean)
			mVT2 <- apply(VT^2, 2, mean)

			# Return the mll #
			(1 + k) / 2 * n * log(2 * pi * sigma2) +
			1 / sigma2 * (
				(1 + k) / 2 * n * VND^2 +
				1 / 2 * sum(CO^2) * sum(VS^2) +
				sum(CO) * VND * sum(VS) -
				t(VS) %*% VT %*% CO -
				n * VND * sum(mVT) +
				n / 2 * sum(mVT2)
			)
		},
		c(rep(0.5, k), 0, sd(VT[, 1]), VT[, 1])
	)

	# Format results #
	x <- list(
		VT = VT,
		coefficients = matrix(rep(c(0, f$estimate[1:k]), n), byrow = TRUE, ncol = 1 + k), 
		VND = matrix(rep(f$estimate[k + 1], 1 + k), ncol = 1 + k),
		VS = matrix(rep(f$estimate[(k + 3):(k + 3 + n - 1)], 1 + k), ncol = 1 + k),
		sigma = matrix(rep(sqrt(f$estimate[k + 2]^2), 1 + k), ncol = 1 + k)
	)
	x
}

occ <- function (VT, method = "reml")
{
	# Get parameters #
	VT <- as.matrix(VT)
	n <- nrow(VT)
	if (n < 2) {
		stop("Less than two rows (ROIs) in the VT matrix")
	}
	k <- ncol(VT) - 1
	if (k < 1) {
		stop("Less than two columns (scans) in the VT matrix")
	}
	if (method != "ols" & method != "ref" & method != "reml") {
		warning(paste("method = \"", method, "\" is not supported. Using \"reml\"", sep = ""))
		method <- "reml"
	}

	# Estimate #
	x <- switch(method,
		ref = occ.ref(VT, n, k),
		ols = occ.ols(VT, n, k),
		reml = occ.reml(VT, n, k)
	)

	# Format results #
	x$fitted.values <- matrix(rep(x$VND, n), byrow = TRUE, ncol = 1 + k) + x$VS * (1 - x$coefficients)
	x$residuals <- VT - x$fitted.values
	colnames(x$coefficients) <- colnames(x$VND) <- colnames(x$VS) <- colnames(x$sigma) <- colnames(x$fitted.values) <- colnames(x$residuals) <- colnames(VT)
	rownames(x$coefficients) <- rownames(x$VS) <- rownames(x$fitted.values) <- rownames(x$residuals) <- rownames(VT)
	rownames(x$VND) <- "VND"
	rownames(x$sigma) <- "sigma"

	x$call <- match.call()
	class(x) <- "occ"
	x
}

plot.occ <- function (x, ...)
{
	n <- nrow(x$VS)
	k <- ncol(x$VS) - 1
	VTmax = max(x$VT, x$fitted.values, na.rm = TRUE)

	# Prepare the plot #
	plotcols = plotrows = ceiling(sqrt(n + 1))
	while (plotrows * (plotcols - 1) >= n + 1) {
		plotcols = plotcols - 1
	}
	par(mfrow = c(plotrows, plotcols))

	# Plot each roi #
	for (roi in 1:n) {
		fitted = rbind(x$VND, x$fitted.values[roi,] - x$VND)
		colnames(fitted) = paste(colnames(fitted), " (", round(x$coefficients[roi,] * 100), "%)", sep="")
		barplot(fitted, space=1, ylim=c(0, 1.1 * VTmax), col = c("gray25", "gray75"), ylab = "Volume of distribution")
		lines(c(0.85, 2.15), rep(x$VT[roi, 1], 2))
		for (scan in 1:k) {
			lines(c(scan * 2 + 0.85, scan * 2 + 2.15), rep(x$VT[roi, scan + 1], 2))
		}
		title(paste(rownames(x$VS)[roi], "occupation"))
	}

	# Finish the plot #
	frame()
	legend("topleft", "Observed total volume of distribution", bty = "n", lty = 1)
	legend("left", c("Free specific volume of distribution", "Non-displaceable volume of distribution"), bty = "n", fill = c("gray75", "gray25"))
	par(mfrow = c(1, 1))
}

print.occ <- function (x, ...)
{
	cat("\nCall:\n")
	print(x$call)
	cat("\nNeuroreceptor occupancy coefficients:\n")
	print(signif(x$coefficients, 6))
	cat("\n")
}

summary.occ <- function (object, ...)
{
	class(object) <- "summary.occ"
	object
}

print.summary.occ <- function (x, ...)
{
	cat("\nCall:\n")
	print(x$call)
	cat("\nResiduals:\n")
	print(signif(summary(c(x$residuals)), 4))
	cat("\nNeuroreceptor occupancy coefficients:\n")
	print(round(x$coefficients, 4))
	cat("\nNon-displaceable volumes of distribution:", signif(x$VND, 4), "\n")
	cat("\nSpecific volumes of distribution:\n")
	print(signif(x$VS, 4))
	cat("\nResidual standard errors:", signif(x$sigma, 4), "\n\n")
}

