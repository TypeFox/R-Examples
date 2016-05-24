toMatrix <-
function(expr.data, gene.sets, flip.gene.sets = FALSE){
  if(is.matrix(gene.sets) | is.array(gene.sets) | is.data.frame(gene.sets) | is.table(gene.sets) ){
    gene.sets <- as.matrix(gene.sets)
    }

  if(is.list(gene.sets)){
    if(flip.gene.sets){
      gene.sets <- flipListStruct(gene.sets)
      gene.sets <- listTOclMatrix(gene.sets)

    }
    else{
      gene.set.matrix <- listTOclMatrix(rownames(expr.data),gene.sets)
      gene.set.size <- dim(gene.set.matrix)
      colnames(gene.set.matrix) <- names(gene.sets)
      rownames(gene.set.matrix) <- rownames(expr.data)
      gene.sets <- gene.set.matrix
    
    }
  }
  return(gene.sets)
}

#############

listTOclMatrix <-
function(gene.vector, GO.list){
	num.classes <- length(GO.list)
	result <- list()
	for(i in 1:num.classes){
		result[[i]] <- gene.vector%in%GO.list[[i]]
	}
	out <- 1*(matrix(unlist(result),ncol=num.classes))
	return(out)
}

#############

flipListStruct <- function(data){
    list_len <- length(data)
    max_val = 0
    ind_vals <- rep(0,list_len)
    loop_val <- 1
    for(k in 1:list_len){
        if(length(data[[k]])> 0){
            if(length(data[[k]])> 1 || data[[k]] > 0) {
                ind_vals[loop_val] <- k
                loop_val <- loop_val + 1
                max_val <- max(max_val, max(data[[k]]) )
            }
        }
        
    }
    ind_vals <- ind_vals[ind_vals > 0]
    out <- list(NULL)
    out[[max_val+1]] = 1
    out[[max_val+1]] = NULL 
    for(k in ind_vals){
        for(l in data[[k]] ){  
            out[[l]] <- c(out[[l]], k)
        }
    }
    out 
}

#############