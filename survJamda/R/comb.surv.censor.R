comb.surv.censor <-
function(geno.files,index,surv.data)
{
	surv = censor = NULL

	for (i in index)
		if (i == 1){
			surv = c(surv,surv.data[[1]][1:nrow(get(geno.files[i]))])
			censor = c(censor,surv.data[[2]][1:nrow(get(geno.files[i]))])
		}
		else{
			curr.ind = 0
			for(j in 2:i)
			      curr.ind = curr.ind+nrow(get(geno.files[j-1]))

			surv = c(surv,surv.data[[1]][(curr.ind+1):(curr.ind+nrow(get(geno.files[i])))])
			censor = c(censor,surv.data[[2]][(curr.ind+1):(curr.ind+nrow(get(geno.files[i])))])
		}		
		
	return(list (surv = surv, censor = censor))
}

