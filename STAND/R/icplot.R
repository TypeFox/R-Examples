icplot <-
function(surv, time = as.numeric(names(surv)), xrange = NA, lines.only = FALSE, 
	XLAB = "Time", YLAB = "Probability", LTY = 1, ...)
{
	k <- length(surv)
	if(length(time) != k)
		stop("length of surv and time must be the\n\tsame")
	if(time[k] == Inf)
		stop("time value = Inf, cannot plot it")
	if(lines.only == FALSE) {
		if(is.na(xrange)) {
			xrange <- range(c(time, 0))
		}
		plot(xrange, c(0, 1), type = "n", xlab = XLAB, ylab = YLAB, ...
			)
	}
	x <- rep(0, 2 * k + 1)
	y <- rep(1, 2 * k + 1)
	for(j in 1:(k - 1)) {
		y[(2 * j + 1):(2 * j + 2)] <- surv[j]
		x[(2 * j):(2 * j + 1)] <- time[j]
	}
	y[2 * k + 1] <- surv[k]
	x[(2 * k):(2 * k + 1)] <- time[k]
	lines(x, y, lty = LTY)
}

