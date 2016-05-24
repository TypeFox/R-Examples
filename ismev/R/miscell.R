# This file contains the following functions:
# identity  q.form  mrl.plot

# "identity"<-
# function(x)
# x

"q.form"<-
function(d, m)
{
#
# ancillary routine
# evaluates quadratic forms
#
	t(as.matrix(d)) %*% m %*% as.matrix(d)
}

"mrl.plot"<-
function(data, umin = min(data), umax = max(data) - 0.1, conf = 0.95, nint = 
	100)
{
#
# function to produce empirical mean residual life plot
# as function of threshold.
# confidence intervals included as well.
#
	x <- xu <- xl <- numeric(nint)
	u <- seq(umin, umax, length = nint)
	for(i in 1:nint) {
		data <- data[data > u[i]]
		x[i] <- mean(data - u[i])
		sdev <- sqrt(var(data))
		n <- length(data)
		xu[i] <- x[i] + (qnorm((1 + conf)/2) * sdev)/sqrt(n)
		xl[i] <- x[i] - (qnorm((1 + conf)/2) * sdev)/sqrt(n)
	}
	plot(u, x, type = "l", xlab = "u", ylab = "Mean Excess", ylim = c(min(
		xl[!is.na(xl)]), max(xu[!is.na(xu)])))
	lines(u[!is.na(xl)], xl[!is.na(xl)], lty = 2)
	lines(u[!is.na(xu)], xu[!is.na(xu)], lty = 2)
}

