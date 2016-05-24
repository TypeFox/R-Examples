
addAIC.fun <-function(mlePP, covariatesAdd, startAdd=NULL,modSim=FALSE,...)
{
	covariates<-mlePP@covariates
	ncovAdd<-dim(covariatesAdd)[2]
	if (is.null(startAdd)) startAdd<-rep(0,ncovAdd)
	nancoef<-names(as.list(mlePP@fullcoef))
	lastcoef<-nancoef[length(nancoef)]
	jj<-match(lastcoef, paste('b', c(0:50), sep=''))#since the first position is 0, j is the following covariate
	newlist<-list(0)
	names(newlist)<-paste('b', jj, sep='')

	AICcurrent<-AIC(mlePP,...)
	AICnew<-NULL
	tind<-mlePP@tind
	if (ncovAdd>=1) 
	{
		for (j in c(1:ncovAdd))
		{
			newlist[[1]]<-startAdd[j]
			nstart<-c(as.list(mlePP@coef),newlist)
		      aux<-update(mlePP, covariates=cbind(covariates, covariatesAdd[,j]), 
			start=nstart,modCI=FALSE,modSim=TRUE, dplot=FALSE)
		      AICnew[j]<-AIC(aux,...)
		}
	

	namcovariates<-dimnames(covariatesAdd)[[2]]
	if (is.null(namcovariates)) namcovariates<- paste('New Covariate',c(1:ncovAdd))
	AICnew<-matrix(AICnew,ncol=1,dimnames=list(namcovariates,'AIC'))

	posminAIC<-which.min(AICnew)

	if (modSim==FALSE)
	{
		 cat(fill=T)
             cat(' Initial model. AIC: ', round(AICcurrent,3), fill=TRUE)
		 cat(' Initial model adding covariates     ', fill=TRUE)
		 print(round(AICnew,3))
             cat(fill=T)
             cat('The best covariate to add is ', namcovariates[posminAIC], fill=TRUE)

	}
	}
	else stop('There are no covariates to add')
	newlist[[1]]<-startAdd[posminAIC]
	return(list(AICadd=AICnew,posminAIC=posminAIC, namecov=namcovariates[posminAIC], AICcurrent=AICcurrent, newCoef=newlist))

}

