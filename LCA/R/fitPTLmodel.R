fitPTLmodel <- function(x,nPairs=10000){
	x <- (x-apply(x,MARGIN=1,min))/(apply(x,MARGIN=1,max)-apply(x,MARGIN=1,min))

	SampleList1 <- sample(1:ncol(x),size=nPairs,replace=TRUE)
	SampleList2 <- sample(1:ncol(x),size=nPairs,replace=TRUE)
	prune <- which(SampleList1==SampleList2)
	nPairs <- nPairs - length(prune)
	SampleList1 <- SampleList1[-prune]
	SampleList2 <- SampleList2[-prune]

	PTL_parameter_tables <- mapply(function(a,b)getPTLparams(x[,a],x[,b]),a=SampleList1,b=SampleList2)

	parameter_dists <- unlist(PTL_parameter_tables[1,])
	parameter_alphas <- unlist(PTL_parameter_tables[3,])
	parameter_betas <- unlist(PTL_parameter_tables[2,])
	parameter_gammas <- unlist(PTL_parameter_tables[4,])

	failed <- which(is.na(parameter_betas))
	cat(paste(sum(!is.na(parameter_betas)),"models successfully fitted out of",nPairs,"\n"))
	if(sum(!is.na(parameter_betas))<2) stop("Rerun model calibration (possibly with higher nPairs), no successful fits this time.")
	parameter_dists <- parameter_dists[-failed]
	parameter_alphas <- parameter_alphas[-failed]
	parameter_betas <- parameter_betas[-failed]
	parameter_gammas <- parameter_gammas[-failed]

	parameter_distOrder <- order(parameter_dists,decreasing=FALSE)
	parameter_dists <- parameter_dists[parameter_distOrder]
	parameter_alphas <- parameter_alphas[parameter_distOrder]
	parameter_betas <- parameter_betas[parameter_distOrder]
	parameter_gammas <- parameter_gammas[parameter_distOrder]

	

	list(alpha=lm(alphas ~ dists, data=data.frame(alphas=parameter_alphas,dists=parameter_dists)),beta=lm(formula = betas ~ dists + I(dists^2) + I(dists^3), data = data.frame(betas=parameter_betas,dists=parameter_dists)),gamma=lm(gammas ~ dists, data=data.frame(gammas=parameter_gammas,dists=parameter_dists)))

}



