plotnorm <- function(
		mean = 0, sd = 1,
		lower = -Inf, upper = Inf,
		n = 750, xlim, xlab = "", ylab = "", 
		groups, lwd = 2, lty = 1,
		main = paste("N(", mean, ",", sd, ")", sep = ""),
		scales = list(y = list(draw = FALSE)),
		col = trellis.par.get('superpose.polygon')$col,
		add.text = TRUE,
		digits = 3,
		cex = 1,
		...)
{
	if (missing(xlim)) {
		xlim <- mean + c(-1, 1) * 4.5 * sd
	}


	x <- seq(xlim[1], xlim[2], by= diff(xlim)/n)
		y <- dnorm(x, mean = mean, sd = sd)

		if (missing(groups)) { groups = ((x > lower) + (x > upper)) }

	middle.prob = round(pnorm(upper, mean, sd) - pnorm(lower, mean, sd), digits)
		upper.prob = round(1 - pnorm(upper, mean, sd), digits)
		lower.prob = round(pnorm(lower, mean, sd), digits)

		xyplot(y ~ x, xlab = xlab, ylab = ylab, lwd = lwd, lty = lty, scales = scales, groups = groups,
				main= main,
				panel = function(x, y, ...) {
				panel.xyplot(x, y, type = 'h' , col = col, ...)
				panel.xyplot(x, y, type = 'l', lwd = 3)
				if (is.finite(upper) && !is.finite(lower) && add.text && FALSE) {
				grid.text(x = .98, y = unit(1, 'npc')-unit(5, 'lines'), 
					bquote(P(X>.(upper)) == .(upper.prob)),
					gp = gpar(cex = cex),
					just = 'right'
					)
				}
				if (is.finite(lower) && add.text && FALSE) {
				grid.text(x = .02, y = unit(1, 'npc')-unit(5, 'lines'), 
					bquote(P(X<.(lower)) == .(lower.prob)), 
					gp = gpar(cex = cex),
					just = 'left'
					)
				}
				if (is.finite(upper) && is.finite(lower) && add.text) {
				grid.text(x = unit(mean(c(lower, upper)), 'native'), y = unit(0, 'native')+unit(1, 'lines'), 
					bquote(.(middle.prob)), gp = gpar(cex = cex),
					just = 'center'
					)
				}
				},
				...)
}

