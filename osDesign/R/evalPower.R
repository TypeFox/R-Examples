evalPower <-
function(waldTest,
									    keep,
									    resNames=NULL)
{	
	## Useful objects
	##
	nDesigns <- dim(waldTest)[2]
	p        <- dim(waldTest)[3]
	
	##
	value <- matrix(NA, nrow=nDesigns, ncol=p, dimnames=resNames)
	##
	for(i in 1:nDesigns)
	{
		##
		tempP  <- waldTest[keep[,i],i,]
		value[i,] <- apply(tempP, 2, mean) * 100
	}
	
  ##
  return(value)
}
