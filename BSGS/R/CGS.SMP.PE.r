CGS.SMP.PE = function(CGS.SMP.Output)
{
	
	r.samples = CGS.SMP.Output$r.samples
	r.est = apply(r.samples, 1, mean)

	beta.samples.rm.0 = CGS.SMP.Output$beta.samples
	beta.samples.rm.0[r.samples == 0] = NA
	beta.est = apply(beta.samples.rm.0, 1, mean, na.rm=T)
	beta.est[is.nan(beta.est)] = 0 

	sigma2.est = mean(CGS.SMP.Output$sigma2.samples)

	list(beta.est = beta.est, r.est = r.est, sigma2.est = sigma2.est)
}




MSE.CGS.SMP = function(Output, Y, X)
{
	pe = CGS.SMP.PE(Output)

	mse = sum( (Y-X%*%cbind(pe$beta.est))^2 )/length(Y)

	return(mse)
}

TCR.TPR.FPR.CGS.SMP = function(Output, True.r, Critical.Point)
{
	pe = CGS.SMP.PE(Output)
	r.est = 1*(pe$r.est>Critical.Point)
	TCR = sum(True.r == r.est)/length(True.r)
	TPR = sum(r.est[True.r==1])/sum(True.r)
	FPR = sum(r.est[True.r==0])/sum(True.r==0)

	list(TCR = TCR, TPR = TPR, FPR = FPR)
}
