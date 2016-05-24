# Automatically generated from all.nw using noweb

findAvailNonInform <- function(ped, avail){

  ## trim persons who are available but not informative b/c not parent
  ## by setting their availability to FALSE, then call findUnavailable()
  ## JPS 3/10/14 add strings check in case of char ids
  pedData <- data.frame(id=ped$id, father=ped$findex, 
                        mother=ped$mindex, avail=avail, stringsAsFactors=FALSE )
  
  checkParent <- is.parent(pedData$id, pedData$father, pedData$mother)
  
  for(i in 1:nrow(pedData)){
    
    if(checkParent[i]==FALSE & avail[i]==TRUE & 
       all(ped$affected[i]==0, na.rm=TRUE)) {

      ## could use ped$affected[i,] if keep matrix
      
        fa <- pedData$id[pedData$father[i]]
        mo <- pedData$id[pedData$mother[i]]
        if(avail[pedData$id==fa] & avail[pedData$id==mo])
          {
            pedData$avail[i] <- FALSE
          }
      }
  }

  idTrim <- findUnavailable(ped, pedData$avail)
  return(idTrim)
} 

