outCSV <- function(block.obj, namesCol = NULL, file.names = NULL, digits = 2, ...){

  ## takes block, assignment, or diagnose object
  if(!is.null(block.obj$blocks)){ 
    block.obj <- block.obj$blocks
  }
  if(!is.null(block.obj$assg)){ 
    block.obj <- block.obj$assg
  }
  
  for(i in 1:length(block.obj)){
    tab <- block.obj[[i]]
    nm <- names(block.obj)[i]
    tab[,ncol(tab)] <- 
    tab[,ncol(tab)] <- round(tab[,ncol(tab)], digits)

    if(!is.null(namesCol)){
      names(tab) <- namesCol
    }
    if(is.null(file.names)){
    	file.name <- paste("Group", nm, ".csv", sep="")
    }else{
    	file.name <- paste(file.names[[i]], ".csv", sep = "")
    	}
    
    write.csv(tab, file = file.name, ...)
  }
}
