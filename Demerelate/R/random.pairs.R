random.pairs <- function(tab.all, allele.column, pairs)
		{
	
     ran.pairs <- lapply(seq(1:2),function(x)
       
    {	
      
      ran.pairsA <- sample(c(tab.all[,allele.column],tab.all[,allele.column+1]),pairs,replace=TRUE)
			ran.pairsB <- sample(c(tab.all[,allele.column],tab.all[,allele.column+1]),pairs,replace=TRUE)
			individual <- as.character(paste(rep(tab.all[1,2],pairs),"random",seq(1:pairs)))
		  population <- rep(as.character(tab.all[1,2]),pairs)
	  	data.frame(individual,population,ran.pairsA,ran.pairsB)
	    
      })
     
			
 return(ran.pairs)

		}

			
