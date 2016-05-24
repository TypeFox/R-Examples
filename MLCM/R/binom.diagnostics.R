binom.diagnostics <- function(obj, nsim = 200, type = "deviance") 
	{
	n <- length(fitted(obj))
	d <- obj$obj$data[, -1]
	ys <- matrix(rbinom(n * nsim, 1, fitted(obj)), 
			ncol = nsim)
	res <- sapply(seq_len(nsim), function(x, obj){
#		ys <- rbinom(n, 1, fitted(obj))
		d <- cbind(resp = ys[, x], d)
		lnk <- obj$obj$family$link
		br <- glm(resp ~ . - 1, binomial(link = lnk), d)
		rs <- residuals(br, type = type)
		rsd <- sort(rs)
		fv.sort <- sort(fitted(br), index.return = TRUE)
		rs <- rs[fv.sort$ix]
		rs <- rs > 0
		runs <- sum(rs[1:(n-1)] != rs[2:n])
		list(resid = as.vector(rsd), NumRuns = runs)
	}, obj = obj
	) 
 	fres <- list(NumRuns = sapply(seq(2, length(res), 2), 
 		function(x) res[[x]]))
 	fres$resid <- t(do.call("cbind", 
 		list(sapply(seq(1, length(res), 2), 
 		function(x) res[[x]]))))
 	fres$resid <- apply(fres$resid, 2, sort)
 	fres$Obs.resid <- residuals(obj$obj, type = type)
	rs <- residuals(obj$obj, type = type)
	fv.sort <- sort(fitted(obj), index.return = TRUE)
	rs <- rs[fv.sort$ix]
	rs <- rs > 0
	obs.runs <- sum(rs[1:(n-1)] != rs[2:n])
	nr <- sum(fres$NumRuns > obs.runs) 
	fres$ObsRuns <- obs.runs
	fres$p <- 1 - nr/nsim
	class(fres) <- c("mlcm.diag", "list")
	fres
}

plot.mlcm.diag <- function(x, alpha = 0.025, 
		breaks = "Sturges", ...) {
	nsim <- dim(x$resid)[1]
	n <- dim(x$resid)[2]
	par(mfrow = c(1, 2))
	plot(sort(x$Obs.resid), (1:n - 0.5)/n, 
		ylab = "Cumulative Density Function",
		xlab = "Deviance Residuals", ...)
	lines(x$resid[alpha * nsim, ], (1:n-0.5)/n, 
			 col = "blue")
	lines(x$resid[(1 - alpha) * nsim, ], (1:n-0.5)/n, 
			 col = "blue")
	hist(x$NumRuns, xlab = "Number of Runs", main = "",
		breaks = breaks)
	abline(v = x$ObsRuns, lwd = 2)
	invisible()
}

