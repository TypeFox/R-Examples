plotCI <-
function(ests, errors, labels = "", zeroline = TRUE, incol = "gray", outcol = "blue", 
	main.title = "Statistical Significance Plot", axes = FALSE, xlab = "", pch = 16, ...) {

	xcorrs = 1:length(ests)
	errors[is.na(errors)] <- 0
	print(errors)
	lowerbounds = ests - 1.96 * errors
	upperbounds = ests + 1.96 * errors
	blackorblue <- function(lower, upper, incol, outcol) {
		if (lower < 0 && 0 < upper) {
			incol
		} else outcol
	}
	col = mapply(blackorblue, lowerbounds, upperbounds, incol, outcol)
	plot(xcorrs, ests, axes = axes, xlab = xlab, col = col, pch = pch, ...)
	axis(2, las = 1)
	box()
	segments(x0 = xcorrs, y0 = lowerbounds, y1 = upperbounds, col = col, lwd = 2, lend = 1)
	
	text(xcorrs, upperbounds, labels, cex = 0.7, col = col, pos = 3, offset = .1)
	
	if (zeroline) 
		abline(h = 0, lty = "dashed")

	title(main.title)

}
