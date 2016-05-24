delManual = function(candidates, observations, nDiff)

{
	plot(candidates)
	title(main=paste("select manually stations to delete (",nDiff,"stations 		to go...)"))
	points(coordinates(observations),pch=19,cex=0.7)
	obs = as.data.frame(coordinates(observations))
	ind = seq(1,length(obs$x))
	obs$ind = ind
	nm = NULL
	# Loop over the 20 stations
	for (i in 1:nDiff)
  	{                                                                    
  	nm = identify(obs, n=1, pos=F, plot=F)
    obs = obs[seq(1,length(obs$x)) != nm,]
  	del = as.data.frame(coordinates(observations))[setdiff(ind,obs$ind),]
  	plot(candidates)
  	title(main=paste("select manually stations to delete (",nDiff-i,"stations to go...)"))
    points(obs,pch=19,cex=0.7)
    points(del,pch="X",col="red",cex=1.2)
  }
	if (i == nDiff) title(sub=paste("stations deleted manually!"),cex.sub = 		1.5, font.sub = 3, col.sub = "red")
  return(observations[-setdiff(ind,obs$ind),])
}

