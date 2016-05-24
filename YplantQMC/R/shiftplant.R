shiftplant <- function(plant, delX, delY, delZ){
  
  n <- plant$nleaves
  
  for(i in 1:n){
    
    plant$leaves[[i]]$XYZ[,1] <- plant$leaves[[i]]$XYZ[,1] + delX
    plant$leaves[[i]]$XYZ[,2] <- plant$leaves[[i]]$XYZ[,2] + delY
    plant$leaves[[i]]$XYZ[,3] <- plant$leaves[[i]]$XYZ[,3] + delZ
  
  }
  
  hasstems <- !all(is.na(plant$stems))
  if(hasstems){
     
    for(j in 1:length(plant$stems)){
      if(!is.null(plant$stems[[j]])){
      plant$stems[[j]]$xyz$from[1] <- plant$stems[[j]]$xyz$from[1] + delX
      plant$stems[[j]]$xyz$from[2] <- plant$stems[[j]]$xyz$from[2] + delY
      plant$stems[[j]]$xyz$from[3] <- plant$stems[[j]]$xyz$from[3] + delZ
      
      plant$stems[[j]]$xyz$to[1] <- plant$stems[[j]]$xyz$to[1] + delX
      plant$stems[[j]]$xyz$to[2] <- plant$stems[[j]]$xyz$to[2] + delY
      plant$stems[[j]]$xyz$to[3] <- plant$stems[[j]]$xyz$to[3] + delZ
      }
    }
    
    for(j in 1:length(plant$branches)){
      if(!is.null(plant$branches[[j]])){
      plant$branches[[j]]$xyz$from[1] <- plant$branches[[j]]$xyz$from[1] + delX
      plant$branches[[j]]$xyz$from[2] <- plant$branches[[j]]$xyz$from[2] + delY
      plant$branches[[j]]$xyz$from[3] <- plant$branches[[j]]$xyz$from[3] + delZ
      
      plant$branches[[j]]$xyz$to[1] <- plant$branches[[j]]$xyz$to[1] + delX
      plant$branches[[j]]$xyz$to[2] <- plant$branches[[j]]$xyz$to[2] + delY
      plant$branches[[j]]$xyz$to[3] <- plant$branches[[j]]$xyz$to[3] + delZ
      }
    }
    
    for(j in 1:length(plant$petioles)){
      if(!is.null(plant$petioles[[j]])){
      plant$petioles[[j]]$xyz$from[1] <- plant$petioles[[j]]$xyz$from[1] + delX
      plant$petioles[[j]]$xyz$from[2] <- plant$petioles[[j]]$xyz$from[2] + delY
      plant$petioles[[j]]$xyz$from[3] <- plant$petioles[[j]]$xyz$from[3] + delZ 
      
      plant$petioles[[j]]$xyz$to[1] <- plant$petioles[[j]]$xyz$to[1] + delX
      plant$petioles[[j]]$xyz$to[2] <- plant$petioles[[j]]$xyz$to[2] + delY
      plant$petioles[[j]]$xyz$to[3] <- plant$petioles[[j]]$xyz$to[3] + delZ 
      }
    }  
    
  }
  
  
  plant$leaftipcoor[,1] <- plant$leaftipcoor[,1] + delX
  plant$leaftipcoor[,2] <- plant$leaftipcoor[,2] + delY
  plant$leaftipcoor[,3] <- plant$leaftipcoor[,3] + delZ
  
  plant$leafbasecoor[,1] <- plant$leafbasecoor[,1] + delX
  plant$leafbasecoor[,2] <- plant$leafbasecoor[,2] + delY
  plant$leafbasecoor[,3] <- plant$leafbasecoor[,3] + delZ
  
return(plant)  
}





