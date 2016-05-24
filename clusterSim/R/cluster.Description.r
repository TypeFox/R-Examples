cluster.Description <-function(x, cl,sdType="sample")
{
modalValue<-function(t)
{		
	len<-length(t)
	ile<-0
	aggregated<-array(NA,c(len,2))
	for(ii in 1:len)
	{
		found<-FALSE
		if (ile!=0)
		for (kk in 1:ile)
		{
			if (as.double(t[ii])==as.double(aggregated[kk,1]))
			{
				found<-TRUE
				aggregated[kk,2]<-aggregated[kk,2]+1
			}
		}
		if (! found)
		{
			ile<-ile+1
			aggregated[ile,1]<-t[ii]	
			aggregated[ile,2]<-1	
		}				
	}	
	if (ile==0)
		NA
	else
	if(ile==1)
		aggregated[1,1]
	else
	{
		aggregated<-aggregated[1:ile,]
		maxCount<-max(aggregated[,2])
		if (sum(aggregated[,2]==maxCount)==1)
			aggregated[aggregated[,2]==maxCount,1]		
		else
			NA				
	}
}
	if (sdType!="sample" && sdType!="population")
		stop("sdType parameter should be one of two values: sample or population")
	if(is.null(dim(x))){
    dim(x)<-c(length(x),1)
	}
	m<-ncol(x)
	n<-length(cl)
	k<-max(cl)
	result<-array(0,c(k,m,5))
	for(i in 1:k)
	{
		for (j in 1:m)
		{
			result[i,j,1]=mean(x[cl==i,j])
			result[i,j,2]=sd(x[cl==i,j])
			if (sdType=="population")
			if (sum(cl==i)==1)
			{
				result[i,j,2]=0				
			}
			else
			{
				result[i,j,2]=sd(x[cl==i,j])*sqrt((sum(cl==i)-1)/sum(cl==i))
			}
			
			result[i,j,3]=median(x[cl==i,j])
			result[i,j,4]=mad(x[cl==i,j])
		}
		if (sum(cl==i)==1)
			result[i,,5]<-x[cl==i,]
		else
		{
			t<-x[cl==i,]
			if (m==1)
			{
				result[i,j,5]<-modalValue(t[j])				
			}
			else
			for (j in 1:m)
			{
				result[i,j,5]<-modalValue(t[,j])
			}
		}
	}
	result
}


