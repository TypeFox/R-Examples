FindGenes<-function(DataLimma,names=NULL){
	
	FoundGenes=list()
	
	if(is.null(names)){
		for(j in 1:length(DataLimma)){
			names[j]=paste("Method",j,sep=" ")	
		}
	}
	
	nrclusters=length(DataLimma[[1]])
	
	
	for(j in 1:nrclusters){
		FoundGenes[[j]]=list()
	}
	
	for(i in 1:length(DataLimma)){ #i == method
		
		for(j in 1:nrclusters){ #j == cluster
			if(!(is.na(DataLimma[[i]][[j]]))[1]){
				
				tempgenes=DataLimma[[i]][[j]]$Genes$TopDE$ID		
				
				
				if(!(is.null(tempgenes))){
					tempgenes=tempgenes[!(is.na(tempgenes))]
					for(k in 1:length(tempgenes)){
						if(!(tempgenes[k] %in% names(FoundGenes[[j]]))){	
							#FoundGenes[[j]][[length(FoundGenes[[j]])+1]]=c()
							FoundGenes[[j]][[length(FoundGenes[[j]])+1]]=names[i]
							names(FoundGenes[[j]])[length(FoundGenes[[j]])]=tempgenes[k]
							
						}
						else if (tempgenes[k] %in% names(FoundGenes[[j]])){
							found=which(names(FoundGenes[[j]])==tempgenes[k])
							FoundGenes[[j]][[found]]=c(FoundGenes[[j]][[found]],names[i])
						}	
					}					
				}	
			}
		}
	}
	
	for(l in 1:nrclusters){
		namesl=names(FoundGenes[[l]])
		if(!(is.null(namesl))){
			
			for(m in l:nrclusters){
				namesm=names(FoundGenes[[m]])	
				Templist=FoundGenes[[m]]
				if(length(namesm)!=0){
					if(l != m){
						
						for(k in 1:length(namesm)){
							
							if (namesm[k] %in% namesl){
								found=which(namesl==namesm[k])
								methods=Templist[[k]]
								del=which(names(FoundGenes[[m]])==namesm[k])
								FoundGenes[[m]][[del]]=c()
								for(a in 1:length(methods)){
									methods[a]=paste(methods[a],"_",m,sep="")
								}
								FoundGenes[[l]][[found]]=c(FoundGenes[[l]][[found]],methods)
							}
							
						}
					}					
				}
			}
		}			
	}
	
	
	for(i in 1:length(FoundGenes)){
		names(FoundGenes)[i]=paste("Cluster",i,sep=" ")
	}
	
	return(FoundGenes)	
}
