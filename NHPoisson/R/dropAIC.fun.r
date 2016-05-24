dropAIC.fun <-function(mlePP,modSim=FALSE,...)
{
	covariates<-mlePP@covariates
	ncov<-dim(covariates)[2]
	AICcurrent<-AIC(mlePP,...)
	AICnew<-NULL
	tind<-mlePP@tind
	if (ncov>(tind==FALSE)) #if there is an intercept ncov can be 1 but otherwise must be at least 2
	{
		for (j in c(1:ncov))
		{
			aux<-update(mlePP, covariates=covariates[,-j], start=as.list(mlePP@coef[-(j+(tind==TRUE))]),modCI=FALSE,modSim=TRUE, dplot=FALSE)
		    AICnew[j]<-AIC(aux,...)
		}

	namcovariates<-dimnames(mlePP@covariates)[[2]]
	if (is.null(namcovariates)) namcovariates<- paste('Covariate',c(1:ncov))
	AICnew<-matrix(AICnew,ncol=1,dimnames=list(namcovariates,'AIC'))

	posminAIC<-which.min(AICnew)

	if (modSim==FALSE)
	{
		 cat(fill=T)
             cat(' Initial model ', round(AICcurrent,3), fill=TRUE)
		 cat(' Initial model deleting covariate ',fill=TRUE)
		 print(round(AICnew,3))
		 cat(fill=T)
             cat('The best covariate to drop is ', namcovariates[posminAIC], fill=TRUE)

	}
	}
	else stop('No test can be carried out since there is only one covariate and no intercept')

	return(list(AICdrop=AICnew,posminAIC=posminAIC,namecov=namcovariates[posminAIC], AICcurrent=AICcurrent))

}













