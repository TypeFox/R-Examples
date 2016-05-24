addManual = function( candidates, observations, nDiff)
{
	observations = as.data.frame(coordinates(observations))
  plot(candidates)
	title(main=paste("locate manually stations to add (",nDiff,"stations to go...)"))
	points(observations,pch=19,cex=0.7)
	nm2=NULL
	for (i in 1:nDiff)
  	{
  	nm2Temp=as.data.frame(locator(n = 1, type="p"))
  	nm2=rbind(nm2,nm2Temp)
  	plot(candidates)
  	title(main=paste("locate manually stations to add (",nDiff-i,"stations to go...)"))
  	points(observations,pch=19,cex=0.7)
  	points(nm2,pch=19,col="green")
  	}
	if (i == nDiff) title(sub=paste("stations added manually!"),cex.sub = 	1.5, font.sub = 3, 	col.sub = "green")
	# Computes the Mukv with the manual network
  add = data.frame(x = c(observations[,1],nm2[,1]), y = c(observations[,2],nm2[,2]))
	coordinates(add) = ~x+y
  return(add)
}
