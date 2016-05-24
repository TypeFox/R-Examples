checkdim<-function(covariates,j)
{
#it checks if a column deletion gives a matrix with only one column since
#in that case the dimnames attribute is lost and must be redefined
	if (dim(covariates)[2]==2)
	{
		auxname<-dimnames(covariates)[[2]][-j]
		covariates<-as.matrix(covariates[,-j])
		dimnames(covariates)<-list(NULL,auxname)
	}
	else
	covariates<-as.matrix(covariates[,-j])

	return(covariates)
}