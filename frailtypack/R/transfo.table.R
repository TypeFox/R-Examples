#
#
#
transfo.table <- function(data,event1,event2,nom,id,recurrentAG){
	
	ni <- as.vector(tapply(data[,id],data[,id],length)) # nb de lignes pour chaque groupe
	
	dataR <- NULL
	dataM <- NULL
	
	for(i in 1:length(ni)){ # pour chaque groupe
		t1cal <- NULL
		if (i==1) ind <- 0
		else ind <- ind + ni[i-1]
		
		tab <- data[(ind+1):(ind+ni[i]),]
		
		if(!recurrentAG){
			if(length(nom)==2){
				k <- 0
				for(j in 1:ni[i]){
					t1 <- k + data[,nom[2]][j+ind]
					k <- t1
					t1cal <- c(t1cal,t1)
				}
			}else{
				k <- 0
				for(j in 1:ni[i]){
					t1 <- k + data[,nom][j+ind]
					k <- t1
					t1cal <- c(t1cal,t1)
				}
			}
			tab <- cbind(tab,t1cal)
		}
		
		last <- tab[nrow(tab),]
		bloc <- tab[which(tab[,event1]==1),]
		bloc2 <- tab[which(tab[,event2]==1),]
		bloc <- rbind(bloc,last)
		bloc2 <- rbind(bloc2,last)
		
		if(!recurrentAG){ # gap-time
			if(length(nom)==2){
				bloc[1,nom[2]] <- bloc[1,ncol(bloc)] 
				bloc2[1,nom[2]] <- bloc2[1,ncol(bloc2)]
				if(nrow(bloc)>1){
					for(j in 2:nrow(bloc)){
						bloc[j,nom[2]] <- bloc[j,ncol(bloc)]-bloc[j-1,ncol(bloc)]
					}
				}
				if(nrow(bloc2)>1){
					for(j in 2:nrow(bloc2)){
						bloc2[j,nom[2]] <- bloc2[j,ncol(bloc2)]-bloc2[j-1,ncol(bloc2)]
					}
				}
			}else{
				bloc[1,nom] <- bloc[1,ncol(bloc)] 
				bloc2[1,nom] <- bloc2[1,ncol(bloc2)]
				if(nrow(bloc)>1){
					for(j in 2:nrow(bloc)){
						bloc[j,nom] <- bloc[j,ncol(bloc)]-bloc[j-1,ncol(bloc)]
					}
				}
				if(nrow(bloc2)>1){
					for(j in 2:nrow(bloc2)){
						bloc2[j,nom] <- bloc2[j,ncol(bloc2)]-bloc2[j-1,ncol(bloc2)]
					}
				}
			}
		}else{ # calendar-time
			bloc[1,nom[1]] <- 0
			bloc2[1,nom[1]] <- 0
			
			if(nrow(bloc)>1){
				for(j in 2:nrow(bloc)){
					bloc[j,nom[1]] <- bloc[j-1,nom[2]]
				}
			}
			if(nrow(bloc2)>1){
				for(j in 2:nrow(bloc2)){
					bloc2[j,nom[1]] <- bloc2[j-1,nom[2]]
				}
			}
		}
		dataR <- rbind(dataR,bloc)
		dataM <- rbind(dataM,bloc2)
	}
	
	return(list(dataR=dataR,dataM=dataM))
	
}
