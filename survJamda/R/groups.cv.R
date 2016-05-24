groups.cv <-
function(n, ngroup, censor)
{
	found = 0
       leave.out <- trunc(n/ngroup)
        groups <- vector("list", ngroup)

	while (!found){
		j = 1
	 	o <- sample(1:n)
        	while (j < ngroup) {
            		jj <- (1 + (j - 1) * leave.out)
            		groups[[j]] <- (o[jj:(jj + leave.out - 1)])
	    		if (sum (censor[groups[[j]]], na.rm=TRUE) ==0)
				break		
	    		else 
				j = j+1	
        	}	
	 	if (j == ngroup){
			groups[[ngroup]] <- o[(1 + (ngroup - 1) * leave.out):n]
			if (sum (censor[groups[[ngroup]]], na.rm=TRUE) > 0)
				break
		}
		
	}
	return(groups)
}

