niche.width <-
function(mat, method = c("shannon","levins")){
    ### Leven 1968 = inverse simpson
    niche.width.levins <- 
    function(mat){ 
       Bi <- c()
       for (i in 1:ncol(mat)){
    	  nij <- mat[,i]
    	  nij <- nij[nij > 0]
    	  pij <- nij/sum(nij)
    	  Bi[i] <- 1/sum(pij^2)
    	  }
    	  Bi <- as.data.frame(t(Bi))
    	  colnames(Bi) <- colnames(mat)
          return(Bi)
    }
    ### Shannon
    niche.width.shannon <- 
    function(mat){ 
       Bi <- c()
       for (i in 1:ncol(mat)){
    	  nij <- mat[,i]
    	  nij <- nij[nij > 0]
    	  pij <- nij/sum(nij)
    	  Bi[i] <- -sum(pij*log(pij))
    	  }
    	  Bi <- as.data.frame(t(Bi))
    	  colnames(Bi) <- colnames(mat)
          return(Bi)
    }
	
    match.arg(method)
    mat <- na.omit(mat)
	if(!is.data.frame(mat)){
	   mat <- as.data.frame(mat)
	}
	if(ncol(mat) < 2){
	   stop("The input data contain less than two species, cannot compute niche overlap")
	}
	match.arg(method)
    if(method == "shannon"){
       result <- niche.width.shannon(mat)
    }
    else {
          if(method == "levins"){
             result <- niche.width.levins(mat)
             }
         }
    return(result)
}

