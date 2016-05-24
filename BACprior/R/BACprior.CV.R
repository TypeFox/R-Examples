BACprior.CV<-function (Y, X, U, omega = c(1, 1.1, 1.3, 1.6, 2, 5, 10, 30, 50, 100, Inf), maxmodels = 150, cutoff = 1e-04, V = 100, criterion = "CVm") 
{
	na.fail(cbind(Y, X, U));
	n = length(X);
	n2 = floor(n/2);
	Criterion = numeric(length(omega));
	sample0 = matrix(0, nrow = V, ncol = (n - n2));
	sample1 = matrix(0, nrow = V, ncol = n2);
	Beta0 = numeric(V);
	betas1 = matrix(0, nrow = V, ncol = length(omega));
	for (v in 1:V) 
	{
		samp1 = sample(1:n, n2);
		samp0 = (1:n)[-samp1];
		sample0[v,] = samp0;
		sample1[v,] = samp1;
		X0 = X[sample0[v,]]
		U0 = U[sample0[v,], ];
		Y0 = Y[sample0[v,]];
		Beta0[v] = BACprior.lm(Y0, X0, U0, omega = Inf, maxmodels = maxmodels, cutoff = cutoff)$results[2];
		X1 = X[sample1[v,]];
		U1 = U[sample1[v,], ];
		Y1 = Y[sample1[v,]];
		betas1[v,] = BACprior.lm(Y1, X1, U1, omega = omega, maxmodels = maxmodels, cutoff = cutoff)$results[, 2];
	}	
	if(criterion == "CVm"){beta0 = mean(Beta0)};
	if(criterion == "CV"){beta0 = Beta0};
	Criterion.Value = colSums((betas1 - beta0)^2/V);
	best = omega[which.min(sqrt(Criterion))];
	plot(omega, Criterion.Value, type = "b", main = "Criterion value according to omega value", 
   	sub = bquote("Best omega = " ~ .(best)), ylab = criterion)
   	abline(v = best, col = "red", lty = 2)
	return(list(best = best, Criterion = Criterion.Value))
}
