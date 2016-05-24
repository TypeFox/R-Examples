"gev" <- 
function(data, block = NA, ...)
{
    n.all <- NA
    if(!is.na(block)) {
        n.all <- length(data)
        if(is.character(block)) {
            times <- as.POSIXlt(attributes(data)$times)
            if(block %in% c("semester", "quarter")) {
                sem <- quart <- times$mon
                sem[sem %in% 0:5] <- quart[quart %in% 0:2] <- 0
                sem[sem %in% 6:11] <- quart[quart %in% 3:5] <- 1
                quart[quart %in% 6:8] <- 2
                quart[quart %in% 9:11] <- 3
            }
            grouping <- switch(block,
                semester = paste(times$year, sem),
                quarter = paste(times$year, quart),
                month = paste(times$year, times$mon),
                year = times$year,
                stop("unknown time period"))
            data <- tapply(data, grouping, max)
        }
        else {
            data <- as.numeric(data)
            nblocks <- (length(data) %/% block) + 1
            grouping <- rep(1:nblocks, rep(block, nblocks))[1:length(data)]
            data <- tapply(data, grouping, max)
        }
    }
    data <- as.numeric(data)
    n <- length(data)
    sigma0 <- sqrt(6 * var(data))/pi
    mu0 <- mean(data) - 0.57722 * sigma0
    xi0 <- 0.1
    theta <- c(xi0, sigma0, mu0)
    negloglik <- function(theta, tmp)
    {
      	y <- 1 + (theta[1] * (tmp - theta[3]))/theta[2]
       	if((theta[2] < 0) || (min(y) < 0))
       	    out <- 1e+06
       	else {
       	    term1 <- length(tmp) * logb(theta[2])
       	    term2 <- sum((1 + 1/theta[1]) * logb(y))
       	    term3 <- sum(y^(-1/theta[1]))
       	    out <- term1 + term2 + term3
       	}
       	out
    }
    fit <- optim(theta, negloglik, hessian = TRUE, ..., tmp = data)
    if(fit$convergence)
        warning("optimization may not have succeeded")
    par.ests <- fit$par
    varcov <- solve(fit$hessian)
    par.ses <- sqrt(diag(varcov))
    out <- list(n.all = n.all, n = n, data = data, block = block, par.ests
       	 = par.ests, par.ses = par.ses, varcov = varcov, converged = 
       	fit$convergence, nllh.final = fit$value)
    names(out$par.ests) <- c("xi", "sigma", "mu")
    names(out$par.ses) <- c("xi", "sigma", "mu")
    class(out) <- "gev"
    out
}

"gumbel" <- 
function(data, block = NA, ...)
{
	n.all <- NA
	data <- as.numeric(data)
        if(!is.na(block)) {
	  n.all <- length(data)
	  if(fg <- n.all %% block) {
              data <- c(data, rep(NA, block - fg))
              warning(paste("final group contains only", fg, "observations"))
          }
          data <- apply(matrix(data, nrow = block), 2, max, na.rm = TRUE)
	}
	n <- length(data)
	sigma0 <- sqrt(6 * var(data))/pi
	mu0 <- mean(data) - 0.57722 * sigma0
	theta <- c(sigma0, mu0)
	negloglik <- function(theta, tmp)
	{
		y <- (tmp - theta[2])/theta[1]
		if(theta[1] < 0)
			out <- 1e+06
		else {
			term1 <- length(tmp) * logb(theta[1])
			term2 <- sum(y)
			term3 <- sum(exp( - y))
			out <- term1 + term2 + term3
		}
		out
	}
        fit <- optim(theta, negloglik, hessian = TRUE, ..., tmp = data)
        if(fit$convergence)
            warning("optimization may not have succeeded")
	par.ests <- fit$par
	varcov <- solve(fit$hessian)
	par.ses <- sqrt(diag(varcov))
	out <- list(n.all = n.all, n = n, data = data, block = block, par.ests
		 = par.ests, par.ses = par.ses, varcov = varcov, converged = 
		fit$convergence, nllh.final = fit$value)
	names(out$par.ests) <- c("sigma", "mu")
	names(out$par.ses) <- c("sigma", "mu")
	class(out) <- "gev"
	out
}

"plot.gev" <- 
function(x, ...)
{
	par.ests <- x$par.ests
	mu <- par.ests["mu"]
	sigma <- par.ests["sigma"]
	if(!("xi" %in% names(par.ests)))
	    xi <- 0
	else xi <- par.ests["xi"]
	if(xi != 0)
	    residuals <- (1 + (xi * (x$data - mu))/sigma)^(-1/xi)
	else residuals <- exp( - (x$data - mu)/sigma)
	choices <- c("Scatterplot of Residuals", "QQplot of Residuals")
	tmenu <- paste("plot:", choices)
	pick <- 1
	while(pick > 0) {
	    pick <- menu(tmenu, title =
                         "\nMake a plot selection (or 0 to exit):")
	    switch(pick,
		   {
		       plot(residuals, ylab = "Residuals",
                            xlab = "Ordering", ...)
		       lines(lowess(1:length(residuals), residuals))
		   },
		   qplot(residuals, ...))
	}
}

"rlevel.gev" <- 
function(out, k.blocks = 20, add = FALSE, ...)
{
	par.ests <- out$par.ests
	mu <- par.ests["mu"]
	sigma <- par.ests["sigma"]
	if(!("xi" %in% names(par.ests)))
	    stop("Use this function after a GEV rather than a Gumbel fit")
	else xi <- par.ests["xi"]
	pp <- 1/k.blocks
	v <- qgev((1 - pp), xi, mu, sigma)
	if(add) abline(h = v)
	data <- out$data
        overallmax <- out$nllh.final
	sigma0 <- sqrt(6 * var(data))/pi
	xi0 <- 0.01
	theta <- c(xi0, sigma0)
	parloglik <- function(theta, tmp, pp, rli)
	{
		mu <- rli + (theta[2] * (1 - ( - logb(1 - pp))^( - theta[
			1])))/theta[1]
		y <- 1 + (theta[1] * (tmp - mu))/theta[2]
		if((theta[2] < 0) | (min(y) < 0))
			out <- 1e+06
		else {
			term1 <- length(tmp) * logb(theta[2])
			term2 <- sum((1 + 1/theta[1]) * logb(y))
			term3 <- sum(y^(-1/theta[1]))
			out <- term1 + term2 + term3
		}
		out
	}
	parmax <- NULL
	rl <- v * c(0.5, 0.6, 0.7, 0.8, 0.85, 0.9, 0.95, 1, 1.1, 1.2,
                    1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3, 3.25, 3.5, 4.5)
	for(i in 1:length(rl)) {
		fit <- optim(theta, parloglik, hessian = FALSE, tmp = data,
                             pp = pp, rli = rl[i])
		parmax <- rbind(parmax, fit$value)
	}
	parmax <-  - parmax
	overallmax <-  - overallmax
	crit <- overallmax - qchisq(0.9999, 1)/2
	cond <- parmax > crit
	rl <- rl[cond]
	parmax <- parmax[cond]
	smth <- spline(rl, parmax, n = 200)
	aalpha <- qchisq(0.95, 1)
	if(!add) {
	    plot(rl, parmax, type = "p", ...)
	    abline(h = overallmax - aalpha/2)
	    abline(v = v)
	    lines(smth)
	}
        ind <- smth$y > overallmax - aalpha/2
	ci <- range(smth$x[ind])
	if(add) {
	    abline(h = ci[1], lty = 2, col = 2)
	    abline(h = ci[2], lty = 2, col = 2)
	}
	as.numeric(c(ci[1], v, ci[2]))
}


