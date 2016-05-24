"getAniA" <- function (A, dataset, ani, anipar) 
{
	if(ani$angle[dataset] != "MA" && !ani$useparperp) {
	   dA <- nrow(A)
	   if(ani$angle[dataset] == "PAR") 
		gamma <- 2
	   if(ani$angle[dataset] == "PERP")
		gamma <- -1
	    for(i in 1:dA) {
		for(j in 1:dA) {
	        if(ani$anifunc[i] == "const" || ani$anifunc[i] == "exp") 
		   A[j,i] <- A[j,i] * (1+(gamma *
		   anipar[ani$parperm[[i]][1]])) 
		if(ani$anifunc[i] == "expvar") 
		   A[j,i] <- A[j,i] * (1 + (gamma * (anipar[ ani$parperm[[i]][1]]
		   - anipar[ ani$parperm[[i]][2]]))) 
			
	     }
	   }
        }
	A
}