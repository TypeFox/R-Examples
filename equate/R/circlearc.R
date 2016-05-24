#----------------------------------------------------------------
# Circle-arc equating functions

#----------------------------------------------------------------
# The main circle-arc function
# type is only included as an argument so it doesn't
# get passed through ... to the linear function

circlearc <- function(x, y, type, method,
	lowp = c(min(scales(x)), min(scales(y))),
	midp = c(mean(x), NA), highp = c(max(scales(x)),
	max(scales(y))), chainmidp = "mean",
	simple = TRUE, verbose = FALSE, ...) {

	if(missing(y))
		y <- margin(x, 2)
	if(margins(y) < margins(x))
		x <- margin(x)

	xscale <- scales(x)
	if(is.na(midp[2])) {
		if(method == "none")
			midp[2] <- mean(y)
		else if(method == "chained") {
			slope2 <- ifelse(chainmidp == "mean", 1,
				sd.freqtab(y)/sd.freqtab(y, 2))
			midp[2] <- mean(y) + slope2*(mean(x, 2) -
				mean(y, 2))
		}
		else
			midp <- linear(x, y, type = "mean", method = method,
				verbose = T, ...)$synth$mean
	}

	xp <- c(lowp[1], midp[1], highp[1])
	yp <- c(lowp[2], midp[2], highp[2])

	slope <- (yp[3] - yp[1])/(xp[3] - xp[1])
	intercept <- yp[1] - slope*xp[1]
	yps <- yp - (xp*slope + intercept)
	if(simple) {
		cent <- c(xcenter = (xp[3]^2 - xp[1]^2)/(2*(xp[3] - xp[1])),
			ycenter = ((xp[1]^2)*(xp[3] - xp[2]) -
				(xp[2]^2 + yps[2]^2)*(xp[3] -
				xp[1]) + (xp[3]^2)*(xp[2] - xp[1]))/
				(2*(yps[2]*(xp[1] - xp[3]))))
		r <- radius(xp[1], 0, cent)
	}
	else {
		cent <- center(xp, yp)
		r <- radius(xp[1], yp[1], cent)
	}

	yx <- circ(xscale, xp, yp, yps, intercept, slope, cent,
		r, simple)

	if(verbose) {
		out <- list(x = x, y = y,
			concordance = data.frame(scale = xscale, yx = yx),
			simple = simple)
		if(method == "chained") out$chainmidp <- chainmidp
		out <- c(out, list(...)[names(list(...)) %in%
			c("lts", "internal", "w")])
		out$coefficients <- c(intercept = intercept,
			slope = slope, cent, r = r)
		out$points <- data.frame(rbind(xp, yp, yps),
			row.names = c("x", "y", "ystar"))
		colnames(out$points) <- c("low", "mid", "high")
	}
	else out <- yx

	return(out)
}

#----------------------------------------------------------------
# Circle-arc conversion of x to y

circ <- function(x, xp, yp, yps, intercept, slope,
	cent, r, simple = TRUE) {

	index <- x >= xp[1] & x <= xp[3]
	out <- slope*x + intercept
	out[index] <- circle(x[index], yps, cent, r) +
		if(simple & yps[2]) out[index] else 0
	names(out) <- NULL

	return(out)
}

#----------------------------------------------------------------
# Circle functions

circle <- function(x, yps, cent, r) {
	
	if(yps[2] < 0)
		out <- cent[2] - sqrt((r^2) -
			(x - cent[1])^2)
	else if(yps[2] > 0)
		out <- cent[2] + sqrt((r^2) -
			(x - cent[1])^2)
	else out <- x
	names(out) <- NULL
	
	return(out)
}

center <- function(x, y) {
	a <- rbind(c(2*(x[1] - x[3]), 2*(y[1] - y[3])),
		c(2*(x[2] - x[3]), 2*(y[2] - y[3])))
	b <- c(x[1]^2 - x[3]^2 + y[1]^2 - y[3]^2,
		x[2]^2 - x[3]^2 + y[2]^2 - y[3]^2)
	out <- solve(a, b)
	names(out) <- c("xcenter", "ycenter")

	return(out)
}

radius <- function(x, y, center) {

	out <- sqrt((center[1] - x[1])^2 +
		(center[2] - y[1])^2)
	names(out) <- NULL
	return(out)
}

#----------------------------------------------------------------
