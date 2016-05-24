LRTpv.fun <-
function(mlePP)
{
	covariates<-mlePP@covariates
	ncov<-dim(covariates)[2]
	LRTpv<-NULL
	tind<-mlePP@tind
	if (ncov>(tind==FALSE)) #if there is an intercept ncov can be 1 but otherwise must be at least 2
	{

		for (j in c(1:ncov))
		{
			llikR<--update(mlePP, covariates=covariates[,-j], start=as.list(mlePP@coef[-(j+(tind==TRUE))]),
			modCI=FALSE,modSim=TRUE, dplot=FALSE)@min
                  difdev <- 2*(-mlePP@min - llikR)
     			LRTpv[j]<-1-pchisq(difdev, 1)
		}
	}
	else stop('No test can be carried out since there is  not enough covariates')

	namcovariates<-dimnames(mlePP@covariates)[[2]]
	if (is.null(namcovariates)) namcovariates<-paste('Covariate',c(1:ncov))
	LRTpv<-matrix(LRTpv,ncol=1,dimnames=list(namcovariates,'p-values'))

      cat(fill=T)
      cat('  The p-values of the LRT comparing the initial model and the model without the covariate ',  fill=TRUE)
	print(round(LRTpv,3))
	return(LRTpv)

}

