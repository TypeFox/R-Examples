# determination of P(k,N,w)
pval <- function(k,N,w)
{	
	
	return((k/w-N-1)*b(k,N,w)+2*Gb(k,N,w))
		
}

# helper function
b<-function(k,N,w)
{
	return(choose(N,k)*w^k*(1-w)^(N-k))	
}

# helper function
Gb<-function(k,N,w)
{	
	sum<-0
	for(i in k:N)
	{	
		sum <- sum + b(i,N,w)	
	}
	
	return(sum)
}

# If two significant overlapping windows were found, these windows are
# merged. If the windows do not overlap, two different windows are stored 
# in a list 
listadapt <- function(lcur,lnew)
{
	if(length(lcur)==0)
	{
		lcur=lnew	
		return(lcur)
	}
	else
	{		
		if(lnew[[1]][1]<=lcur[[length(lcur)]][2])
		{
			lcur[[length(lcur)]][2]<-lnew[[1]][2]
			if(lcur[[length(lcur)]][3]>lnew[[1]][3])
			{
				lcur[[length(lcur)]][3] <- lnew[[1]][3]
			}
			return(lcur)
		}	
		else
		{
			lcur<-append(lcur,lnew)
			return(lcur)
		}
			
		
	}
	
}
# This method searches for data accumulations by shifting a window with
# window size <w> across the data and deciding at each position if there 
# is a data accumulation. To test this, a scan statistic with significance
# level <sign.level> is used. 
scanStatistic <- function(vect, w=0.25, sign.level=0.1)
{	
	temp<-vect
	vect <-unlist(vect)
	vsort <- sort(vect)
	N <- length(vect)
	range <- (max(vect)) - (min(vect))
	windowsize <- range*w	
	N <- length(vect)
	binarizeddata<-temp
	res<-list()
	lcur<-list()
	
	# shift a fixed window over the data
	# the window is moved from point to point
	for(i in seq_along(vect))
	{	
		start <- vsort[i]
		stop <- vsort[i] + windowsize
		
		k <- length(vect[(vect >= start) & (vect <= stop)])
		
		p <- pval(k,N,w)
		
		if(p>1)
		{
			p=0.99	
		}
		
		if(p<=sign.level & p>0 & k >= (N*w-1) & k > 2)
		{
			res <- listadapt(res,list(c(start,stop,p)))
		}
		
	}
	
	
	# if no accumulation for a fixed <sign.level> was found, the 
	# binarization is rejected, and we search for a accumulation
	# with a higher sign.level.
	if(length(res)==0)
	{ 
		while(TRUE)
		{
			sign.level=sign.level+0.05
			
			if(sign.level>2)
			{
				binarizeddata<-(sapply(vect,function(x) 0))  
				return(list(bindata=binarizeddata,thresholds=NA,reject=TRUE))
    			
			}
		
			for(i in seq_along(vect))
			{
				start <- vsort[i]
				stop <- vsort[i] + windowsize
		
				k <- length(vect[(vect >= start) & (vect <= stop)])
		
				p <- pval(k,N,w)
				
				if(p>1)
				{
					p=0.99	
				}
				
				if(p<=sign.level & p>0 & k >= (N*w-1) & k > 2)
				{	
					#res <- append(res,list(c(start=start,stop=stop,pval=p)))
					res <- listadapt(res,list(c(start,stop,p)))
				}
		
			}
			if(length(res)!=0)
				break
				
				
			
		}
		reject<-TRUE	
		
	}
	else
	{
		reject<-FALSE
	}
	
	
	# search the window with the smallest sign.level.
	# this window is used for the binarization
	min=1000
	ind=0
	for(i in seq_along(res))
	{
		if(res[[i]][3]<min)
		{
			ind=i	
			min=res[[i]][3]
		}
	}
	
	# are more points on the left or on the right side
	# of the window? Based on this, the binarization is performed
	bigger <- length(vect[vect > res[[ind]][2]]) 
	smaller <- length(vect[vect < res[[ind]][1]])
	
	if(bigger > smaller)
	{
		threshold<-res[[ind]][2]
		
		small<-tail(vsort[vsort<=threshold],n=1)
		big<-vsort[vsort>threshold][1]
		
		thres<-(big+small)/2
		
		for(i in seq_along(vect))
		{
			if(vect[i]<=threshold)
			{
				binarizeddata[i]<-0
			}
			else
			{
				binarizeddata[i]<-1
			}
		}	
			
	}
	else
	{
		threshold<-res[[ind]][1]
		
		small<-tail(vsort[vsort<threshold],n=1)
		big<-vsort[vsort>=threshold][1]		
		
		thres<-(big+small)/2
		
		for(i in seq_along(vect))
		{
			if(vect[i]>=threshold)
			{
				binarizeddata[i]<-1
			}
			else
			{
				binarizeddata[i]<-0
			}
		}
		
	}
	
	return(list(bindata=binarizeddata,thresholds=as.numeric(thres),reject=reject))
	
}
