
# 
# Author: robin
###############################################################################

.FDRS_Selection<-function(S,q,n)
{
	X=sort(S,decreasing = TRUE,index.return=TRUE)
	som=0;
	#for (i in 1:length(X$x)){som=som+1/i}
	p=q/(2*length(X$x))*(1:length(X$x))
	t=qnorm(p,lower.tail=FALSE);
	last=1;
	for (i in 1:length(S))
	{   
		if (X$x[i]-t[i]>0)
		{
			last=i
		}
	}
	#print(last)
	return(X$ix[1:last])
}



.FDR_Selection<-function(S,q)
{
    X=sort(S,decreasing = TRUE,index.return=TRUE)
	som=0;
	#for (i in 1:length(X$x)){som=som+1/i}
	p=q/(2*length(X$x))*(1:length(X$x))
	t=qnorm(p,lower.tail=FALSE);
	last=1;
	for (i in 1:length(S))
	{	
		if (X$x[i]-t[i]>0)
		{
			last=i
		}
	}
	#print(last)
	return(X$ix[1:last])
}

.FAN_Selection<-function(T,Cmoins,n)
{
	X=sort(T,decreasing = TRUE,index.return=TRUE)
	stat=array(0,length(T))
	for (i in 1:length(T))
	{	
		stat[i]=min(Cmoins[X$ix[1:i]])*(sum(X$x[1:i]^2))^2/(4*i/n+sum(X$x[1:i]^2))
	}
	#print(last)
	return(X$ix[1:which.max(stat)])
}

