aa.conservation_first <-
function(modal, p, seqlength){
	conserve <- 0
	for(i in 1:seqlength){
		if(modal[i] >= p){
			conserve <- conserve + 1
		}
	}
	return(conserve)
}
