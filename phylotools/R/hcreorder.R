#### Function hreorder as part of R package phylotools
#### By Jinlong Zhang  <Jinlongzhang01@gmail.com>
#### Institute of Botany, the Chinese Academy of Sciences, Beijing ,China
#### Nov- 01-2010

hcreorder <- 
function(x)
{
    height = x$height
    merge  = x$merge
    
    ind.height = order(height)
    height2 = height[ind.height]
    merge2 = merge[ind.height,]
    
	## Whether rows of matrix have value greater than limit
    greater <- 
	function(x,limit){
        res00 <- rep(FALSE, nrow(x))
        for(i in 1:nrow(x)){
            if(any(x[i,] > limit)) {
	    	    res00[i] = TRUE
	    	}
        }
	    res00
	}
	
	## Find a matched rownumber
    findrow <- 
    function(xxx, aaa){
        res <- rep(0, nrow(xxx))
        for(i in 1:nrow(xxx)){
    		if(all(xxx[i,] == aaa)){
    			res[i] <- 1
    		}
    	}
    which(res > 0)
    }
	
	merge4 <- merge2
    for(i in which(greater(merge2, 0))){
			## rownumber in matrix merge 
			ind <- which((height2[i] == height)&(merge2[i, 1] == merge[,1])&(merge2[i, 2] == merge[, 2]))
			## rows that contain value > 0 
			selecte1 <- merge[ind,]
			## elements < -1 will not be changed 
			if(selecte1[1] > 0){
			    merge4[i,1] <- findrow(merge2, merge[selecte1[1],])
			}
			if(selecte1[2] > 0){
			    merge4[i,2] <- findrow(merge2, merge[selecte1[2],])
			}
    }
	
    x$height <- height2
    x$merge  <- merge4
    return(x)
}

