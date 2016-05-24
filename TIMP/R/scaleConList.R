"scaleConList" <- 
function (conList, specList) 
{
	maxDim <- max(unlist(lapply(conList, ncol)))
	listProd <-  vector("list", maxDim)
	for(i in 1:length(conList)) {
	    for(j in 1:maxDim) {	
		if(j <= ncol(conList[[i]])) {
		  dtemp <-  abs(conList[[i]][,j] %*% t(specList[[i]][,j]))
		  listProd[[j]] <- append(listProd[[j]], dtemp)
	        }
            } 
        }
	compMax <- unlist(lapply(listProd, max))
	maxall <- max(compMax)
	factScale <- sqrt(compMax / maxall) 
	for(i in 1:length(conList)) {
		conListMax <- apply(conList[[i]], 2, max)
		conList[[i]] <- t(t(conList[[i]]) * (factScale / conListMax))
	}
	conList
}

