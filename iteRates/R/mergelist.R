mergelist <-
function(pG){
	startcombs<-list(pG[[1]])

	startrows<-length(startcombs[[1]][,1])
	startcols<-length(startcombs[[1]][1,])
	addcombs<-pG[[2]]
	multiplier<-dim(addcombs[[1]])
	
	newcomb<-array(dim=c(startrows*multiplier[1],startcols+multiplier[2]))

	k<-1
		for (i in 1:startrows){
		for (j in 1:multiplier[1]){
			newcomb[k,]<-c(startcombs[[1]][i,],addcombs[[i]][j,])
			k<-k+1
						}
		}
	return(newcomb)	
	}

