lhc_constructCoEffDataSet <-
function(LHCRESULTFILE,PARAMNAME,PARAMETERS)
{
	coEffData<-NULL
	coEffHeaders<-NULL

	for(m in 1:length(PARAMETERS))
	{	
		if(PARAMETERS[m]!=PARAMNAME)
		{
			coEffData<-cbind(coEffData,LHCRESULTFILE[,PARAMETERS[m]])
			coEffHeaders<-cbind(coEffHeaders,PARAMETERS[m])
		}
		colnames(coEffData)<-c(coEffHeaders)
	}
	return(coEffData)
}

