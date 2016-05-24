#
# function for calculating which components of a mixture model should be merged
#
# - takes a model from Mclust and returns a list of components to merge
#
##################################################################################

mergeComponents <- function(model,overlap){

	mergelist <- list()

	if(model$G > 2){ # otherwise there is no need to do any merging
	
		mergelist[[1]] <- 1
		mergelist[[2]] <- model$G
		
		table <- c()
		for(class in 1:model$G){
			if(sum(model$classification==class)==0){
				table <- rbind(table,rep(0,model$G))
			}
			if(sum(model$classification==class)==1){
				table <- rbind(table,model$z[model$classification==class,])
			}
			if(sum(model$classification==class)>1){
				table <- rbind(table,colMeans(model$z[model$classification==class,]))
			}
		}
	
		truthArray <- array(sapply(table,FUN=function(x){x > overlap}),dim=dim(table))
		
		if(min(truthArray[c(1,2),c(2,1)]) == 1){
			if((sum(truthArray[c(1,2),model$G]) + sum(truthArray[model$G,c(1,2)])) == 0){
				mergelist[[1]] <- c(1,2)
				if(model$G > 3){
					if(min(truthArray[c(1,2,3),c(3,2,1)]) == 1){
			                        if((sum(truthArray[c(1,2,3),model$G]) + sum(truthArray[model$G,c(1,2,3)])) == 0){
	                		                mergelist[[1]] <- c(1,2,3)
							if(model$G > 4){
								if(min(truthArray[c(1,2,3,4),c(4,3,2,1)]) == 1){
						                        if((sum(truthArray[c(1,2,3,4),model$G]) + sum(truthArray[model$G,c(1,2,3,4)])) == 0){
						                                mergelist[[1]] <- c(1,2,3,4)
		                        				}
						                }
							}
			                        }
			                }
				}
			}
		}
		G <- model$G
		if(min(truthArray[c(G,G-1),c(G-1,G)]) == 1){
                        if((sum(truthArray[c(G,G-1),1]) + sum(truthArray[1,c(G,G-1)])) == 0){
                                mergelist[[2]] <- c(G,G-1)
                                if(model$G > 3){
                                        if(min(truthArray[c(G,G-1,G-2),c(G-2,G-1,G)]) == 1){
                                                if((sum(truthArray[c(G,G-1,G-2),1]) + sum(truthArray[1,c(G,G-1,G-2)])) == 0){
                                                        mergelist[[2]] <- c(G,G-1,G-2)
                                                        if(model$G > 4){
                                                                if(min(truthArray[c(G,G-1,G-2,G-3),c(G-3,G-2,G-1,G)]) == 1){
                                                                        if((sum(truthArray[c(G,G-1,G-2,G-3),1]) + sum(truthArray[1,c(G,G-1,G-2,G-3)])) == 0){
                                                                                mergelist[[2]] <- c(G,G-1,G-2,G-3)
                                                                        }
                                                                }
                                                        }
                                                }
                                        }
                                }
                        }
                }

	}
	else{
		mergelist[[1]] <- 1
		mergelist[[2]] <- model$G
	}
	mergelist
}

