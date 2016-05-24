BACprior.boot<-function (Y, X, U, omega = c(1, 1.1, 1.3, 1.6, 2, 5, 10, 30, 50, 100, Inf), maxmodels = 150, cutoff = 1e-04, B = 100) 
{
	na.fail(cbind(Y, X, U));
	MSE = numeric(length(omega));
	beta0 = BACprior.lm(Y, X, U, omega = Inf, maxmodels = maxmodels, cutoff = cutoff)$results[2];
	dat = cbind(Y,X,U);
	betas1f = function(dat, i)
	{
		Y = dat[i,1];
		X = dat[i,2];
		U = dat[i, -c(1,2)];
		BACprior.lm(Y, X, U, omega = omega, maxmodels = maxmodels, cutoff = cutoff)$results[, 2];
	}
	betas1 = boot(data = dat, statistic = betas1f, R = B)$t;
	MSE = colSums((betas1 - beta0)^2/B);
	best = omega[which.min(sqrt(MSE))];
	plot(omega, sqrt(MSE), type = "b", main = "Estimated sqrt-MSE according to omega value", 
   	sub = bquote("Best omega = " ~ .(best)))
   	abline(v = best, col = "red", lty = 2)
	return(list(best = best, MSE = MSE))
}