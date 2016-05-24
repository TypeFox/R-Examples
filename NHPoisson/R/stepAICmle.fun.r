
stepAICmle.fun<-function(ImlePP, covariatesAdd=NULL, startAdd=NULL, direction='forward',...)
{

if ((direction!='backward')&(direction!='forward')&(direction!='both'))
	stop('Direction argument not valid, it must be one  of: backward, forward,
or both')


ncov<-dim(ImlePP@covariates)[2]
if ( is.null(ncov)) ncov<-0
if ( (is.null(dimnames(ImlePP@covariates)[[2]]))&(ncov>=1) ) dimnames(ImlePP@covariates)<-list(NULL, paste('Covariate',c(1:ncov)) )


i<-1
dropnames<-NULL
dropaddnames<-NULL
addnames<-NULL

tind<-ImlePP@tind
if(direction=='backward')
{

	AICold<-Inf
	AICnew<--Inf
	mlePPaux<-ImlePP
	covariates<-ImlePP@covariates

	while ((AICold>=AICnew)&( (dim(covariates)[2]+(tind==TRUE)) >1))
	{
		cat(fill=T)
		cat(fill=T)
		cat('Step Backward', i, fill=T)
		cat(fill=T)
		auxdrop<-dropAIC.fun(mlePPaux,...)
		AICold<-auxdrop$AICcurrent
		j<-auxdrop$posminAIC
		AICnew<-auxdrop$AICdrop[j]

		if (AICold>=AICnew)
		{
			covariates<-checkdim(covariates,j=j)
			start<-as.list(mlePPaux@coef[-j])
			names(start)<-paste('b',(c(0:(length(start)-1))+(tind==FALSE)),sep='')
			mlePPaux<-update(mlePPaux, covariates=covariates, start=start,modCI=FALSE,modSim=TRUE, dplot=FALSE)
			cat('This drop improves the model', fill=T)
			dropnames<-c(dropnames,auxdrop$namecov)
		}
		else cat('This drop  does not improve the model', fill=T)
      	i<-i+1

	} #end while
	cat(fill=T)
	cat(fill=T)
	cat('Final model', fill=T)
	cat(fill=T)
	cat('The covariates dropped from the initial model are: ', fill=T)
	print(dropnames)
	cat(fill=T)
}

if(direction=='forward')
{

	if (is.null(covariatesAdd)) stop('For the forward direction, a matrix   must be provided in the argument covariatesAdd')
	covariatesAdd<-as.matrix(covariatesAdd)

	AICold<-Inf
	AICnew<--Inf
	mlePPaux<-ImlePP
	ncovAdd<-dim(covariatesAdd)[2]
	if (is.null(dimnames(covariatesAdd)[[2]])) dimnames(covariatesAdd)<- list(NULL,paste('New Covariate',c(1:ncovAdd)) )
	if (is.null(startAdd)) startAdd<-rep(0,ncovAdd)

	while ((AICold>=AICnew)&(ncovAdd >0))
	{
		cat(fill=T)
		cat(fill=T)
		cat('Step Forward', i, fill=T)

		auxadd<-addAIC.fun(mlePPaux, covariatesAdd=covariatesAdd,startAdd=startAdd,...)
		AICold<-auxadd$AICcurrent
		j<-auxadd$posminAIC
		AICnew<-auxadd$AICadd[j]

		if (AICold>=AICnew)
		{
			covariates<-cbind(mlePPaux@covariates,covariatesAdd[,j])
			addnames<-c(addnames,auxadd$namecov)
			dimnames(covariates)<-list(NULL, c(dimnames(mlePPaux@covariates)[[2]],auxadd$namecov) )
			start<-c(as.list(mlePPaux@coef), auxadd$newCoef)
			mlePPaux<-update(mlePPaux, covariates=covariates, start=start,modCI=FALSE,modSim=TRUE, dplot=FALSE)

			cat('This covariate improves the model', fill=T)
			covariatesAdd<-checkdim(covariatesAdd,j=j)
			ncovAdd<-dim(covariatesAdd)[2]
			startAdd<-startAdd[-j]
			
		}
		else cat('This covariate does not improve the model', fill=T)
		i<-i+1

	} #end while
	cat(fill=T)
	cat(fill=T)
	cat('Final model', fill=T)
	cat(fill=T)
	cat('The covariates added to the initial model are: ', fill=T)
	print(addnames)
	cat(fill=T)


}

if(direction=='both')
{
	
	if (is.null(covariatesAdd)) stop('For the both direction, a matrix   must be provided in the argument covariatesAdd')
	covariatesAdd<-as.matrix(covariatesAdd)
	mlePPaux<-ImlePP
	covariates<-ImlePP@covariates
	changeF<-1
	changeB<-1
	mlePPaux<-ImlePP
	ncovAdd<-dim(covariatesAdd)[2]
	if (is.null(dimnames(covariatesAdd)[[2]])) dimnames(covariatesAdd)<- list(NULL, paste('New Covariate',c(1:ncovAdd)) )
	if (is.null(startAdd)) startAdd<-rep(0,ncovAdd)

	while (max(changeF, changeB)>0)  
	{
# Forward step
	if (ncovAdd >0)
	{
		cat(fill=T)
		cat(fill=T)
		cat('Step Forward', i, fill=T)
		auxadd<-addAIC.fun(mlePPaux, covariatesAdd=covariatesAdd,startAdd=startAdd,...)
		AICold<-auxadd$AICcurrent
		j<-auxadd$posminAIC
		AICnew<-auxadd$AICadd[j]

		if (AICold>=AICnew)
		{
			covariates<-cbind(mlePPaux@covariates,covariatesAdd[,j])
			dropaddnames<-c(dropaddnames,auxadd$namecov)
			names(dropaddnames)[length(dropaddnames)]<-'added'
			dimnames(covariates)<-list(NULL, c(dimnames(mlePPaux@covariates)[[2]],auxadd$namecov) )
			start<-c(as.list(mlePPaux@coef), auxadd$newCoef)
			mlePPaux<-update(mlePPaux, covariates=covariates, start=start,modCI=FALSE,modSim=TRUE, dplot=FALSE)
			cat('This covariate improves the model', fill=T)
			covariatesAdd<-checkdim(covariatesAdd,j=j)
			startAdd<-startAdd[-j]
			ncovAdd<-dim(covariatesAdd)[2]

		}
		else
		{
			cat('This covariate does not improve the model', fill=T)
			changeF<-0
		}
	}

#Backward step

	if ((dim(covariates)[2]+(tind==TRUE)) >1)
	{
		cat(fill=T)
		cat('Step Backward', i, fill=T)
		auxdrop<-dropAIC.fun(mlePPaux,...)
		AICold<-auxdrop$AICcurrent
		j<-auxdrop$posminAIC
		AICnew<-auxdrop$AICdrop[j]

		if (AICold>=AICnew)
		{
			covariatesAdd<-cbind(covariatesAdd, covariates[,j])
			dropaddnames<-c(dropaddnames,auxdrop$namecov)
			names(dropaddnames)[length(dropaddnames)]<-'dropped'
			dimnames(covariatesAdd)[[2]][dim(covariatesAdd)[2]]<-auxdrop$namecov
			startAdd<-c(startAdd,as.numeric(mlePPaux@coef[j]))
			covariates<-checkdim(covariates,j=j)
			start<-as.list(mlePPaux@coef[-j])
			names(start)<-paste('b',(c(0:(length(start)-1))+(tind==FALSE)),sep='')
			mlePPaux<-update(mlePPaux, covariates=covariates, start=start,modCI=FALSE,modSim=TRUE, dplot=FALSE)
			cat('This drop improves the model', fill=T)
		}
		else
		{
			 cat('This drop  does not improve the model', fill=T)
		       changeB<-0
		}
	}
	i<-i+1
	} #end while
	cat(fill=T)
	cat(fill=T)
	cat('Final model', fill=T)
	cat(fill=T)
	cat('The covariates added and dropped to the initial model are: ', fill=T)
	print(dropaddnames)
	cat(fill=T)


} #end if both

cat(fill=T)
cat('Summary of the model: ', fill=T)
print(summary(mlePPaux))
cat(fill=T)
 return(mlePPaux)
}





