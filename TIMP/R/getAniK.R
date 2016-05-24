"getAniK" <- function (k, dataset, ani, anipar)
{
	if(ani$angle[dataset] != "MA") {
	    for(i in 1:length(k)) {
	        if(ani$anifunc[i] == "exp") 
			k[i] <- k[i] + anipar[ ani$parperm[[i]][2] ]
		if(ani$anifunc[i] == "expvar") 
		        k[i] <- k[i] + anipar[ ani$parperm[[i]][3] ]
			
	    }
	}
	k 
}
