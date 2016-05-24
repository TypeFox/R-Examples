add.col <-
function(inputA, inputB, add, according){
    if(!is.data.frame(inputA)){
	    inputA <- as.data.frame(inputA)
	  }
	  
    if(ncol(inputA)==1){
        inputA <- data.frame(seq(1:nrow(inputA)),inputA)
        colnames(inputA) <- c("col0", colnames(inputA))
     }
	 
    colaccordA <- which(colnames(inputA)==according)
    if (is.null(colaccordA)){
	    stop(paste("There is no column called \"",according,"\"in inputA"))
	}
		
    inputA <- inputA[order(as.character(inputA[,colaccordA])),]
    speciesA <- table(as.character(inputA[,colaccordA]))
    inputA.species.uniq <- names(speciesA)
   
   
    if(!is.data.frame(inputB)){
	    inputB <- as.data.frame(inputB)
	}
	
    if(ncol(inputB)==1){
        inputB <- data.frame(seq(1:nrow(inputA)),inputA)
        colnames(inputB) <- c("col0", colnames(inputB))
    }
    colaccordB <- which(colnames(inputB)==according)
	
    if (is.null(colaccordB)){
	    stop(paste("There is no column called \"",according,"\"in inputB"))
	}
    
    inputB.species <- as.character(inputB[,colaccordB])
    coladd <- which(colnames(inputB)==add)
    
	if (!any(colnames(inputB)==add)){
	    stop(paste("There is no colum called \"", add ,"\"in inputB, please check!"))
	}
    inputB.add <- as.character(inputB[,coladd])
    result.add <- rep(NA, length(inputA.species.uniq))
    
    for (i in 1:length(inputA.species.uniq)){
        for (j in 1:length(inputB.species)){
            if (inputA.species.uniq[i] == inputB.species[j]) 
            result.add[i] <- inputB.add[j];
        }
     }
    
	addeddata <- rep(result.add, speciesA)
    result <- data.frame(inputA, addeddata)
    colnames(result) <- c(colnames(inputA), add)
   
    if(any(is.na(add))){
        rownum <- is.na(result$add)
        cat("Warning: NA are found in the results, please check!\n")
        error <- result[rownum,]
        error1 <- head(error)
        print(error1)
    }
return(result)
}
