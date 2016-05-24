spreadanalysis <-
function(g, timedelays, numsamples, normalizebyname=FALSE)
{
	numreached <- matrix(NA, nrow=numsamples, ncol=length(timedelays))

	vtos <- sample(length(V(g)), numsamples)

  vnames <- V(g)[vtos]$Name
  
	for (i in 1:numsamples)
	{	  
    thisv <- V(g)[ vtos[i] ]  
    
    whichv <- V(g)[subcomponent(graph=g, v=thisv, mode="out")]
    
		for (j in 1:length(timedelays))
		{
			starttime <- min(whichv$Time)
			
			numreached[i,j] <- length(unique(
			  whichv[whichv$Time < (starttime + timedelays[j])]$Name
          ))	
		}
		
		print(i/numsamples)
	}
	
  if (normalizebyname)
  {
    tempresult <- data.frame(vnames, numreached / length(unique(V(g)$Name)))
  }
  else
  {
    tempresult <- data.frame(vnames, numreached/length(V(g)))    
  }
	names(tempresult) <- c("startvertex",paste(timedelays))
  
	return(tempresult)
}

