#----------------------------------------------------------------
# Linear equating functions

#----------------------------------------------------------------
# The main function

linear <- function(x, y, type = "linear", method = "none",
	ws = -1, wax, way, wbx, wby, lowp = c(min(scales(x)),
	min(scales(y))), midp = c(median(c(lowp[1], highp[1])),
	median(c(lowp[2], highp[2]))), highp = c(max(scales(x)),
	max(scales(y))), sx = highp[1] - lowp[1], sy = highp[2] -
	lowp[2], cx = midp[1], cy = midp[2], internal = TRUE,
	lts = FALSE, verbose = FALSE, ...) {

	if (missing(y))
		y <- margin.table(x, 2)
	if (margins(y) < margins(x))
		x <- margin.table(x, 1)
	xscale <- scales(x)
	
	if (method == "chained") {
		if (type == "mean")
			slope1 <- slope2 <- 1
		else {
			slope1 <- sd.freqtab(y)/sd.freqtab(y, 2)
			slope2 <- sd.freqtab(x, 2)/sd.freqtab(x)
		}
		intercept <- slope1 * (mean(x, 2) -
			slope2 * mean(x) - mean(y, 2)) + mean(y)
		slope <- slope1*slope2
	}
	else {
		if (method == "none") {
			sigmax <- sd.freqtab(x)
			sigmay <- sd.freqtab(y)
			mux <- mean(x)
			muy <- mean(y)
		}
		else {
			synth <- synthetic(x, y, ws, method, internal, lts)
			stats <- synth$synthstats
			if (!lts) {
				sigmax <- stats$sd[1]
				sigmay <- stats$sd[2]
				mux <- stats$m[1]
				muy <- stats$m[2]
			}
			else {
				sigmax <- synth$gamma[1]
				sigmay <- synth$gamma[2]
				mux <- mean(x)
				muy <- mean(y) + synth$gamma[2] *
					(mean(x, 2) - mean(y, 2))
			}
		}
		if (type %in% c("identity", "general linear")) {
			if (missing(wax)) wax <- 0
			if (missing(way)) way <- 0
			if (missing(wbx)) wbx <- 0
			if (missing(wby)) wby <- 0
		}
		else if (type == "mean") {
			if (missing(wax)) wax <- 0
			if (missing(way)) way <- 0
			if (missing(wbx)) wbx <- 1
			if (missing(wby)) wby <- 1
		}
		else if (type == "linear") {
			if (missing(wax)) wax <- 1
			if (missing(way)) way <- 1
			if (missing(wbx)) wbx <- 1
			if (missing(wby)) wby <- 1
		}
		if (cx == "cy") cx <- muy
		if (cy == "cx") cy <- mux
		if (sx == "sy") sx <- sigmay
		if (sy == "sx") sy <- sigmax
		# Put sigma on the scale of s
		# sigmayr <- sigmay*sx/sigmax
		# Then, sigmaxr <- sx
		slope <- (way*sigmay + (1 - way)*sy)/
			(wax*sigmax + (1 - wax)*sx)
		if (is.na(slope)) slope <- 1
		intercept <- (wby*muy + (1 - wby)*cy) -
			slope*(wbx*mux + (1 - wbx)*cx)
	}
	yx <- lin(xscale, intercept, slope)

	if (verbose) {
		out <- list(x = x, y = y, concordance = data.frame(scale =
			xscale, yx = yx), internal = internal, lts = lts,
			coefficients = c(intercept = intercept, slope = slope),
			points = data.frame(low = lowp, mid = midp,
				high = highp, row.names = c("x", "y")))
		if (method != "chained") {
			out$coefficients[c("cx", "cy", "sx", "sy")] <-
				c(cx, cy, sx, sy)
			out$weights <- c(wax = wax, way = way,
				wbx = wbx, wby = wby)
		}
		if (!method %in% c("none", "chained") & !lts)
			out <- c(out, synth)
	}
	else out <- yx

	return(out)
}

#----------------------------------------------------------------
# Linear conversion of x to y

lin <- function(x, intercept, slope) {

	out <- x * slope + intercept
	names(out) <- NULL
	
	return(out)
}

#----------------------------------------------------------------
