ldr.slices <-
function(y, nslices=3)
{	
	endpoints = function(n, nslices)
	{
		# This function is intended to determine the end-points of the slices.
		increment <- floor(n/nslices); 
		if (nslices ==1) return(n);

		if (nslices > 1) 
		{
			ends<- seq(1:nslices)*increment;
			rest <- n%%nslices;
			if (rest==0) return(ends);

			if (rest>0)
			{
				for (i in 1:rest) ends[i]<-ends[i]+i; 
				for (i in (rest+1):nslices) ends[i]<- ends[i]+rest;
				return(ends)
			}
		}
	}
	n <- length(y); indicators <- vector(length=n);
	sorty <- sort(y); bins.y <- vector("list", nslices);
	ends <- endpoints(n, nslices);

	if (nslices==1)
	{
		bins.y[[1]] <- sorty[1:ends[1]]; 
		indicators[1:ends[1]]<-1;
		return(list(bins=bins.y, nslices=nslices, slice.size=n, slice.indicator=indicators[rank(y)]))
	}
	else if (nslices==2)
	{
		bins.y[[1]] <- sorty[1:ends[1]]; 
		bins.y[[2]] <- sorty[(ends[1]+1):ends[2]];
 		indicators[1:ends[1]]<-1;
 		indicators[(ends[1]+1):ends[2]]<-2;

		return(list(bins=bins.y, nslices=nslices, slice.size=diff(c(0,ends)), slice.indicator=indicators[rank(y)]))
	}
	else
	{		
		bins.y[[1]] <- sorty[1:(ends[1]-1)]; 
		indicators[1:(ends[1]-1)]<-1;

		for (i in 2:(nslices-1))
		{ 
			bins.y[[i]] <- sorty[ends[i-1]:(ends[i]-1)];
			indicators[ends[i-1]:(ends[i]-1)] <- i;
		}
		bins.y[[nslices]] <- sorty[ends[nslices-1]:(ends[nslices])];
		indicators[ends[nslices-1]:(ends[nslices])] <- nslices;

		return(list(bins=bins.y, nslices=nslices, slice.size=diff(c(0,ends)), slice.indicator=indicators[rank(y)]))
	}
}
