"loess.demo" <-
function(x, y, span = 2/3, degree = 1, nearest = FALSE, xlim =
	numeric(0), ylim = numeric(0), verbose = FALSE)
{
	# function to demonstrate the locally weighted regression function loess
	# written by Dr. Greg Snow
	# Brigham Young University, Department of Statistics
	# gls@byu.edu now greg.snow@imail.org
	# Modified by Henrik Aa. Nielsen, IMM, DTU (han@imm.dtu.dk)
	miss.xy <- is.na(x) | is.na(y)
	x <- x[!miss.xy]
	y <- y[!miss.xy]
	y <- y[order(x)]
	x <- x[order(x)]
	fit.d <- loess(y ~ x, degree = degree, span = span,
		family = "gaussian", control = loess.control(
		surface = "direct"))
	fit.i <- loess(y ~ x, degree = degree, span = span,
		family = "gaussian")
	est <- list(x = seq(min(x), max(x), len = 500))
	est$y <- predict(fit.i, newdata = data.frame(x = est$
		x))
	xl <- range(x, est$x, xlim)
	xl <- xl + c(-1, 1) * 0.03 * diff(xl)
	yl <- range(y, est$y, fitted(fit.d), ylim)
	yl <- yl + c(-1, 1) * 0.05 * diff(yl)
	fitPlot <- function(x, y, est, fit.d, xl, yl)
	{
		plot(x, y, pch = 3, xlim = xl, ylim = yl)
		lines(x, fitted(fit.d), col = 'red')
		mtext("Exact estimate with linear interpolation between x-values",
			col = 'red', adj = 0.5, line = 0.5)
		lines(est, col = 'blue')
		mtext("Estimate obtained using the default interpolation scheme",
			col = 'blue', adj = 0.5, line = 2)

		NULL
	}
	fitPlot(x, y, est, fit.d, xl, yl)
	repeat {
		x0 <- locator(1)$x
		if(length(x0) < 1)
			break
		if(nearest)
			x0 <- unique(x[abs(x - x0) == min(
				abs(x - x0))])
		if(verbose){
			cat("x0 =", x0, "\n")
                        flush.console()
                      }
		if(span < 1) {
			q <- as.integer(span * length(x))
			d <- sort(abs(x - x0))[q]
		}
		else {
			d <- max(abs(x - x0)) * sqrt(span)
		}
		w <- rep(0, length(x))
		s <- abs(x - x0) <= d
		w[s] <- (1 - (abs(x[s] - x0)/d)^3)^3
		fitPlot(x, y, est, fit.d, xl, yl)
		symbols(x, y, circles = sqrt(w), inches = 0.3,
			add = T, col = 'lightgrey')
		if(degree > 0)
			lines(x, fitted(lm(y ~ poly(x, degree
				), weights = w)), col = 'purple',
				err = -1)
		else {
			##lines(x, fitted(lm(y ~ 1, weights = w)), col = 8, err = -1)
			abline(a = sum(w * y)/sum(w), b = 0,
				col = 'purple')
		}
		abline(v = x0, col = 'green')
		if(x0 - d > xl[1])
			abline(v = x0 - d, col = 'green', lty = 2)
		if(x0 + d < xl[2])
			abline(v = x0 + d, col = 'green', lty = 2)
	}
}

