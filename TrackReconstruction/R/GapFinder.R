GapFinder <-
function(rawdata,timediff=1, timeformat="%d-%b-%Y %H:%M:%S")
{
	dims=dim(rawdata)
	if('DateTime' %in% colnames(rawdata))
	{
		DateTime=strptime(DateTime,format=timeformat,tz="GMT")
		DateTime2=strptime(DateTime2,format=timeformat,tz="GMT")
	}else{
		DateTime=paste(rawdata[1:(dims[1]-1),1],rawdata[1:(dims[1]-1),2])
		DateTime=strptime(DateTime,format=timeformat,tz="GMT")
		DateTime2=paste(rawdata[2:dims[1],1],rawdata[2:dims[1],2])
		DateTime2=strptime(DateTime2,format=timeformat,tz="GMT")
	}
	
	Gaps=matrix(0,1000,4)
	Gaps=as.data.frame(Gaps)
	colnames(Gaps)<-c("Line","TimeDiffinSec","Time1","Time2")
	Counter=1
	
	DIFF=as.double(abs(difftime(DateTime,DateTime2, tz="GMT",units="sec")))
	for(i in 1:(dims[1]-1))
	{
		if(DIFF[i]>timediff)
		{
			Gaps[Counter,1]=i;
			Gaps[Counter,2]=DIFF[i];
			Gaps[Counter,3]=as.character(DateTime[i]);
			Gaps[Counter,4]=as.character(DateTime2[i]);
			Counter=Counter+1;
		}
	}
	Gaps<-Gaps[1:Counter,]
	return(Gaps=Gaps)
}
