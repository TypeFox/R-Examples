plotnetworkslices <-
function(slices, timedeltas)
{
	timestrings <- paste(timedeltas[,1]," - ", timedeltas[,2])
	
	par(mfrow=c(ceiling(sqrt(length(slices))),ceiling(sqrt(length(slices)))))
	par(mar=c(1,1,1,1))
	for (i in 1:length(slices)) 
	{ 
		plottanet(slices[[i]])
		text(0,0,timestrings[i]) 
	}

}

