plot.bigRR <-
function(x, alpha = .36, ...) {
	par(mfrow = c(3, 1))
	plot(x$u, main = 'Shrinkage Estimates', ylab = 'Effect', xlab = 'Index', 
		 col = rgb(1, 0, 0, alpha), type = 'h')
	plot(abs(x$u), main = paste('Intra-class Correlation:', round(x$lambda/(x$lambda + x$phi), digits = 6)), 
	     ylab = '|Effect|', xlab = 'Index', col = rgb(1, 0, 1, alpha), type = 'h')
 	plot(x$leverage, main = 'Leverages', ylab = 'Size', xlab = 'Index', col = rgb(0, 0, 1, alpha), pch = 19, cex = .36)
}

