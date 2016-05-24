niche.overlap <-
function(mat, method = c("levins","schoener","petraitis","pianka","czech","morisita") ){
    match.arg(method)
    mat <- na.omit(mat)
	result <- rep(NA, ncol(mat)*ncol(mat))
    dim(result) <- c(ncol(mat), ncol(mat))
	rownames(result) <- colnames(mat)
	colnames(result) <- colnames(mat)
    for(i in 1:(ncol(mat)-1)){
       for(j in (i+1):ncol(mat)){
    	  if(method == "levins"){
		    result[j,i] <- niche.overlap.pair(mat[,i],mat[,j], method = "levins")
		  }
    	    else {
		        if(method == "schoener"){
     		      result[j,i] <- niche.overlap.pair(mat[,i],mat[,j], method = "schoener")
		          } 
			     else{
    	             if(method == "petraitis") {
		              result[j,i] <- niche.overlap.pair(mat[,i],mat[,j], method = "petraitis")
                      }
		              else{
    	                   if(method == "pianka"){
		                      result[j,i] <- niche.overlap.pair(mat[,i],mat[,j], method = "pianka")
                            }
                            else{		  
    	                        if(method == "czech"){       
		                            result[j,i] <- niche.overlap.pair(mat[,i],mat[,j], method = "czech")
							      }
                                  else{		  
    	                            if(method == "morisita"){    
		                                result[j,i] <- niche.overlap.pair(mat[,i],mat[,j], method = "morisita")
                                       }							  
    	                            }
		                
       						}
		                   
				      }
		        
				 }
		    
			}
        
		}
	
	}
    as.dist(result)
}

