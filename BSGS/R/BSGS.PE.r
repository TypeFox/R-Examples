BSGS.PE = function(BSGS.Output)
{
	
	r.samples = BSGS.Output$r.samples
	r.est = apply(r.samples, 2, mean)

	beta.samples.rm.0 = BSGS.Output$beta.samples
	beta.samples.rm.0[r.samples == 0] = NA
	beta.est = apply(beta.samples.rm.0, 1, mean, na.rm=T)
	beta.est[is.nan(beta.est)] = 0 

	eta.samples = BSGS.Output$eta.samples
	eta.est = apply(BSGS.Output$eta.samples, 2, mean)

	sigma2.est = mean(BSGS.Output$sigma2.samples)

	list(beta.est = beta.est, eta.est = eta.est, r.est = r.est, sigma2.est = sigma2.est)
}


MSE.BSGS = function(Output, Y, X)
{
	pe = BSGS.PE(Output)

	mse = sum( (Y-X%*%cbind(pe$beta.est))^2 )/length(Y)

	return(mse)
}

TCR.TPR.FPR.BSGS = function(Output, True.r, Critical.Point)
{
	pe = BSGS.PE(Output)
	r.est = 1*(pe$r.est>Critical.Point)
	TCR = sum(True.r == r.est)/length(True.r)
	TPR = sum(r.est[True.r==1])/sum(True.r)
	FPR = sum(r.est[True.r==0])/sum(True.r==0)

	list(TCR = TCR, TPR = TPR, FPR = FPR)
}
