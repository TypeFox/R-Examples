det.set.ind <- function(geno.files,train,i)
{
           curr.ind = 0
           ind = NULL
          for (m in i){
	  if(train)	cat (geno.files[m], " ")
		if (m == 1)
                        ind = c(ind,1:nrow(get(geno.files[1])))
                else{
			 curr.ind = 0
                        curr.ind = curr.ind+nrow(get(geno.files[m-1]))

                        ind = c(ind,(curr.ind+1):(curr.ind+nrow(get(
geno.files[m]))))
                }    
	}
     return(ind)
}
