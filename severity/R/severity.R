severity = function(mu0, xbar, sigma, n, alpha)
{
	require(graphics)
	### check inputs ###
	mu0 = as.integer(mu0) # make sure it is an integer
  n = as.integer(n) # make sure it is an integer
  sigma = as.numeric(sigma)
  alpha = as.numeric(alpha)
  if(class(xbar) != "numeric")
  {
    xbar = as.numeric(xbar)
  }
  ### begin severity calculations ###
	r = length(xbar)
	gamma = seq(from = -0.5, to = 1.5, by = 0.05) # discrepancies of interest
	l = length(gamma)
	mu1 = rep(x = mu0, times = l) + gamma
	c_alpha = qnorm(alpha, lower.tail = FALSE) # cut-off (for rejection region)
	sigma_x = sigma / sqrt(n)
	sigma_x_inv = 1 / sigma_x
	d_x0 = sigma_x_inv * (xbar - rep(x = mu0, times = r)) # test statistic
	# delta1 = [sqrt(n) * (mu1 - mu0)] / sigma
	delta1 = sigma_x_inv * gamma # non-centrality parameter
	p = pnorm(d_x0, lower.tail = FALSE) # p-value
	power = pnorm(q = rep(x = c_alpha, times = l) - delta1, lower.tail = FALSE) # power curve
	# "accept" indicates whether to 'accept H0' (accept = 1) OR to 'reject H0' (accept = 0)
	accept = rep(0, times = r) # reject null hypothesis (default)
	for(k in 1:r)
	{
		if(d_x0[k] < c_alpha)
		{
			accept[k] = 1 # accept null hypothesis
		}
	}
	x1 = xbar[which(accept == 0)]
	x2 = xbar[which(accept == 1)]
	m1 = length(x1)
	m2 = length(x2)	
	q1 = matrix(nrow = l, ncol = m1)
	q2 = matrix(nrow = l, ncol = m2)
	sev_rejectH0 = matrix(nrow = l, ncol = m1)
	sev_acceptH0 = matrix(nrow = l, ncol = m2)
	x_reject = matrix(nrow = l, ncol = m1) # for plotting
	x_accept = matrix(nrow = l, ncol = m2) # for plotting
	if((m1 > 0) && (m2 > 0))
	{
		par(mfrow = c(1, 2))
	}
	### reject H0 ###
	if(m1 > 0)
	{
		for(i in 1:m1)
		{
			x_reject[, i] = mu1 # for plotting
			q1[, i] = sigma_x_inv * (rep(x = x1[i], times = l) - mu1)
			sev_rejectH0[, i] = pnorm(q = q1[, i]) # severity
		}
		# plot
		plot(x = x_reject[, 1], y = sev_rejectH0[, 1], type = "l", lty = 1, col = "blue", main = expression(paste("severity curves for inference: ", mu > mu[1], "")), xlab = expression(paste("values of ", mu[1], "")), ylab = "severity / power", xaxt = "n", yaxt = "n", cex.main = 0.90)
		axis(side = 1, at = seq(from = mu1[1], to = mu1[l], by = 0.10), labels = seq(from = mu1[1], to = mu1[l], by = 0.10)) # x-axis
		axis(side = 2, at = seq(from = 0, to = 1, by = 0.05), labels = seq(from = 0, to = 1, by = 0.05)) # y-axis
		lines(x = mu1, y = power, col = "red", lty = 2) # power curve
		if(m1 == 1)
		{
		  if(m2 > 0)
		  {
		    legend(x = mu1[floor(l/2) - 3], y = 0.525, col = c("blue", "red"), lty = c(1, 2), legend = c(as.expression(bquote(bar(x) == .(x1[1]))), "power"), cex = 0.65)
		  }
      else
      {
        legend(x = mu1[l - 10], y = 0.525, col = c("blue", "red"), lty = c(1, 2), legend = c(as.expression(bquote(bar(x) == .(x1[1]))), "power"), cex = 0.65)
      }
		}
    
	}
	if(m1 == 2)
	{
		lines(x = x_reject[, 2], y = sev_rejectH0[, 2], col = "blue", lty = 3)
		if(m2 > 0)
		{
		  legend(x = mu1[floor(l/2) - 3], y = 0.525, col = c("blue", "blue", "red"), lty = c(1, 3, 2), legend = c(as.expression(bquote(bar(x) == .(x1[1]))), as.expression(bquote(bar(x) == .(x1[2]))), "power"), cex = 0.65)
		}
    else
    {
      legend(x = mu1[l - 10], y = 0.525, col = c("blue", "blue", "red"), lty = c(1, 3, 2), legend = c(as.expression(bquote(bar(x) == .(x1[1]))), as.expression(bquote(bar(x) == .(x1[2]))), "power"), cex = 0.65)
    }
	}
	else if(m1 >= 3)
	{
		lines(x = x_reject[, 2], y = sev_rejectH0[, 2], col = "blue", lty = 3)
		lines(x = x_reject[, 3], y = sev_rejectH0[, 3], col = "blue", lty = 4)
		if(m2 > 0)
		{
			legend(x = mu1[floor(l/2) - 3], y = 0.5, col = c("blue", "blue", "blue", "red"), lty = c(1, 3, 4, 2), legend = c(as.expression(bquote(bar(x) == .(x1[1]))), as.expression(bquote(bar(x) == .(x1[2]))), as.expression(bquote(bar(x) == .(x1[3]))), "power"), cex = 0.6)
		}
		else
		{
			legend(x = mu1[l - 10], y = 0.5, col = c("blue", "blue", "blue", "red"), lty = c(1, 3, 4, 2), legend = c(as.expression(bquote(bar(x) == .(x1[1]))), as.expression(bquote(bar(x) == .(x1[2]))), as.expression(bquote(bar(x) == .(x1[3]))), "power"), cex = 0.6)
		}
	}
	### accept H0 ###
	if(m2 > 0)
	{
		for(j in 1:m2)
		{
			x_accept[, j] = mu1 # for plotting
			q2[, j] = sigma_x_inv * (rep(x = x2[j], times = l) - mu1)
			sev_acceptH0[, j] = pnorm(q = q2[, j], lower.tail = FALSE) # severity
		}
		# plot
		plot(x = x_accept[, 1], y = sev_acceptH0[, 1], type = "l", lty = 1, col = "blue", main = expression(paste("severity curves for inference: ", mu <= mu[1], "")), xlab = expression(paste("values of ", mu[1], "")), ylab = "severity / power", xaxt = "n", yaxt = "n", cex.main = 0.90)
		axis(side = 1, at = seq(from = mu1[1], to = mu1[l], by = 0.10), labels = seq(from = mu1[1], to = mu1[l], by = 0.10)) # x-axis
		axis(side = 2, at = seq(from = 0, to = 1, by = 0.05), labels = seq(from = 0, to = 1, by = 0.05)) # y-axis
		lines(x = mu1, y = power, col = "red", lty = 2) # power curve
		if(m2 == 1)
		{
			legend(x = mu1[floor(l/2) - 2], y = 0.275, col = c("blue", "red"), lty = c(1, 2), legend = c(as.expression(bquote(bar(x) == .(x2[1]))), "power"), cex = 0.65)
		}
	}
	if(m2 == 2)
	{
		lines(x = x_accept[, 2], y = sev_acceptH0[, 2], col = "blue", lty = 3)
		legend(x = mu1[floor(l/2) - 2], y = 0.275, col = c("blue", "blue", "red"), lty = c(1, 3, 2), legend = c(as.expression(bquote(bar(x) == .(x2[1]))), as.expression(bquote(bar(x) == .(x2[2]))), "power"), cex = 0.65)
	}
	else if(m2 >= 3)
	{
		lines(x = x_accept[, 2], y = sev_acceptH0[, 2], col = "blue", lty = 3)
		lines(x = x_accept[, 3], y = sev_acceptH0[, 3], col = "blue", lty = 4)
		legend(x = mu1[floor(l/2) - 2], y = 0.25, col = c("blue", "blue", "blue", "red"), lty = c(1, 3, 4, 2), legend = c(as.expression(bquote(bar(x) == .(x2[1]))), as.expression(bquote(bar(x) == .(x2[2]))), as.expression(bquote(bar(x) == .(x2[3]))), "power"), cex = 0.6)
	}
	output = data.frame(sev_rejectH0, sev_acceptH0, power, gamma)
	return(list(accept = accept, p = p, "severity_acceptH0" = sev_acceptH0, "severity_rejectH0" = sev_rejectH0, power = power, "discrepancy" = gamma))
}
