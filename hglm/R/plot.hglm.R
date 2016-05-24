`plot.hglm` <-
	function(x, pch = "+", pcol = 'slateblue', lcol = 2, device = NULL, name = NULL, ...) {
	
	residuals <- x$resid
	fitted.values <- x$fv
	disp.residuals <- x$disp.resid
	disp.fitted.values <- x$disp.fv
	hatvalues <- x$hv
	deviances <- x$dev
	x$nRand <- cumsum(x$RandC)
	#p <- x$p
	#cook.dis <- deviances/(p*sum(deviances))*hatvalues/(1 - hatvalues)**2
	if (is.null(nrow(x$SummVC1))) idx <- c(1, 3, 4) else idx <- 1:4
	for (i in idx) {
		if (i > 1) {
			if (is.null(device)) dev.new() 
		}
		if (i == 1) {
			if (!is.null(device)) pdf(paste(name, i, '.pdf', sep = ''), width = 10, height = 10)
			par(mfrow = c(2, 2), pty = "s", ...)
			loess.fit <- loess.smooth(fitted.values, residuals)
			plot(fitted.values, residuals, xlab = "Fitted Values", 
				 ylab = "Studentized Residuals", pch = pch, col = pcol, bty = "n", main = "Mean Model (a)")
			lines(loess.fit$x, loess.fit$y, col = lcol)
			loess.fit <- loess.smooth(fitted.values, abs(residuals))
			plot(fitted.values, abs(residuals), xlab = "Fitted Values", 
				 ylab = "|Studentized Residuals|", pch = pch, col = pcol, bty = "n", main = "Mean Model (b)")
			lines(loess.fit$x, loess.fit$y, col = lcol)
			qqnorm(residuals, col = pcol, pch = pch, bty = "n", 
				   xlab = "Normal Quantiles", ylab = "Residual Quantiles", main = "Mean Model (c)")
			qqline(residuals, col = lcol)
			hist(residuals, density = 15, xlab = "Studentized Residuals", main = "Mean Model (d)", col = pcol)
			if (!is.null(device)) dev.off()
		}
		else {
			if (i == 2) {
				if (!is.null(device)) pdf(paste(name, i, '.pdf', sep = ''), width = 10, height = 10)
				par(mfrow = c(2, 2), pty = "s", ...)
				loess.fit <- loess.smooth(disp.fitted.values, disp.residuals)
				plot(disp.fitted.values, disp.residuals, xlab = "Fitted Values", 
					 ylab = "Standardized Deviance Residuals", pch = pch, col = pcol, bty = "n", main = "Dispersion Model (a)")
				lines(loess.fit$x, loess.fit$y, col = lcol)
				loess.fit <- loess.smooth(disp.fitted.values, abs(disp.residuals))
				plot(disp.fitted.values, abs(disp.residuals), xlab = "Fitted Values", 
					 ylab = "|Standardized Deviance Residuals|", pch = pch, col = pcol, bty = "n", main = "Dispersion Model (b)")
				lines(loess.fit$x, loess.fit$y, col = lcol)
				qqnorm(disp.residuals, col = pcol, pch = pch, bty = "n", 
					   xlab = "Normal Quantiles", ylab = "Residual Quantiles", main = "Dispersion Model (c)")
				qqline(disp.residuals, col = lcol)
				hist(disp.residuals, density = 15, xlab = "Standardized Deviance Residuals", main = "Dispersion Model (d)", col = pcol)
				if (!is.null(device)) dev.off()
			}
			else {
				if (i == 3) {
					if (!is.null(device)) pdf(paste(name, i, '.pdf', sep = ''), width = 10, height = 5)
					par(mfrow = c(1, 2), pty = "s", ...)
					plot(hatvalues, ylab = "Hat-values", main = "Hat-values", pch = pch, col = pcol, bty = "n")
					plot(deviances, ylab = "Deviances", main = "Deviances", pch = pch, col = pcol, bty = "n")
					if (!is.null(device)) dev.off()
				} else {
					if (!is.null(device)) pdf(paste(name, i, '.pdf', sep = ''), width = 10, height = (length(x$RandC) + 1)*5)
					par(mfrow = c(length(x$RandC) + 1, 2))
					devid <- 1:(length(deviances) - max(x$nRand))
					beta <- var(deviances[devid])/mean(deviances[devid])
					alpha <- mean(deviances[devid])/beta
  					if (length(deviances[devid]) < 5001) max.L <- 10000
  					if (length(deviances[devid]) > 5000) max.L <- length(deviances[devid])*10
  					xx <- rgamma(max.L, alpha, 1/beta) 
  					steps <- floor(max.L/length(deviances[devid]))
  					vec.indx <- steps*(1:length(deviances[devid])) - round(steps/2)
  					x.alt <- sort(xx)[vec.indx]
  					sy <- sort(deviances[devid])
  					plot(x.alt, sy, main = "Mean Model Deviances", col = pcol, pch = pch, bty = "n", ylab = "Deviance Quantiles", xlab="Gamma Quantiles")
					abline(0, 1, col = lcol)
					hist(deviances[devid], density = 15, xlab = "Deviances", main = "Mean Model Deviances", col = pcol)
					devid <- (length(deviances) - max(x$nRand) + 1):(length(deviances) - max(x$nRand) + x$nRand[1])
					beta <- var(deviances[devid])/mean(deviances[devid])
					alpha <- mean(deviances[devid])/beta
  					if (length(deviances[devid]) < 5001) max.L <- 10000
  					if (length(deviances[devid]) > 5000) max.L <- length(deviances[devid])*10
  					xx <- rgamma(max.L, alpha, 1/beta) 
  					steps <- floor(max.L/length(deviances[devid]))
  					vec.indx <- steps*(1:length(deviances[devid])) - round(steps/2)
  					x.alt <- sort(xx)[vec.indx]
  					sy <- sort(deviances[devid])
  					plot(x.alt, sy, main = paste(names(x$SummVC2)[1], "Deviances"), col = pcol, pch = pch, bty = "n", ylab = "Deviance Quantiles", xlab="Gamma Quantiles")
					abline(0, 1, col = lcol)
					hist(deviances[devid], density = 15, xlab = "Deviances", main = paste(names(x$SummVC2)[1], "Deviances"), col = pcol)
					if (length(x$RandC) > 1) {
						for (J in 2:length(x$RandC)) {
							devid <- (length(deviances) - max(x$nRand) + x$nRand[J - 1] + 1):(length(deviances) - max(x$nRand) + x$nRand[J])
							beta <- var(deviances[devid])/mean(deviances[devid])
							alpha <- mean(deviances[devid])/beta
  							if (length(deviances[devid]) < 5001) max.L <- 10000
  							if (length(deviances[devid]) > 5000) max.L <- length(deviances[devid])*10
  							xx <- rgamma(max.L, alpha, 1/beta) 
  							steps <- floor(max.L/length(deviances[devid]))
  							vec.indx <- steps*(1:length(deviances[devid])) - round(steps/2)
  							x.alt <- sort(xx)[vec.indx]
  							sy <- sort(deviances[devid])
  							plot(x.alt, sy, main = paste(names(x$SummVC2)[J], "Deviances"), col = pcol, pch = pch, bty = "n", ylab = "Deviance Quantiles", xlab="Gamma Quantiles")
							abline(0, 1, col = lcol)
							hist(deviances[devid], density = 15, xlab = "Deviances", main = paste(names(x$SummVC2)[J], "Deviances"), col = pcol)
						}
					}
					if (!is.null(device)) dev.off()
				}
			}
		}
	}
}

