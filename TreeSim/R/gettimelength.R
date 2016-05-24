gettimelength <-
function(numbs,n){
		timeperiods <- which(numbs[,3]==n)
		timelength<-vector()
		if (length(timeperiods)>0){
		for (i in 1:length(timeperiods)){
			timelength <- c(timelength, numbs[timeperiods[i],2]-numbs[timeperiods[i]+1,2])
		}	
		}
		timelength
	}

