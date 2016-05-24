percenta <-
function(x,first,last)
{
{
	percentages<-x
	sum<-as.vector(apply(x[,first:last],1,sum))
	for(i in 1:nrow(x))
		{
		percentages[i,]<-x[i,]*100/sum[i]
		} 
}
return(percentages)
}

