`calcpos` <-
function(cross, step, fixpos)
{
  n.chr <- nchr(cross)
  npop <- nind(cross)
  n.mrk <- nmar(cross)
  pos <- list()
  map <- pull.map(cross)

  if (fixpos==0)
      output <- calc.genoprob(cross, step=step)
 
  if ((step==0)&(fixpos>0))
  {
  	pos <- list()
 	for (i in 1:n.chr)
	{
	  pos[[i]] <- rep(map[[i]][1:(n.mrk[i]-1)], each=fixpos+1)+as.vector(sapply(diff(map[[i]]), function(x) return(seq(0, x*fixpos/(fixpos+1), x/(fixpos+1)))))
	  pos[[i]] <- c(pos[[i]], map[[i]][length(map[[i]])])
	  names(pos[[i]]) <- paste("loc", round(pos[[i]], 1), sep="")
	  names(pos[[i]])[seq(1, length(pos[[i]]), fixpos+1)] <- names(map[[i]])
	}
	class(pos) <- "map"
	output <- calc.genoprob2(cross, pos)
  } 
  return(output)
}

