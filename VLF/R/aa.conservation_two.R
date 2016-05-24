aa.conservation_two <-
function(modal1, modal2, p, seqlength){
	both <- modal1 + modal2
	conserve <- 0
	for(i in 1:seqlength){
		if(both[i] >= p){
			conserve <- conserve + 1
		}
	}
	return(conserve)
}
