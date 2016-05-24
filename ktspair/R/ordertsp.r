
##Function used to obtain the order in which the pair of genes should appear (from the best one to the worst one)

ordertsp <- function(delta, gamma){
	order <- order(delta, decreasing=TRUE)
	delta2 <- delta[order]
	gamma2 <- gamma[order]
	unique <- unique(delta2)
	for(i in 1:length(unique)){
		if(length(which(delta2==unique[i]))>1){
			sens <- order(gamma2[which(delta2==unique[i])], decreasing=TRUE)
			order[which(delta2==unique[i])] <- order[which(delta2==unique[i])][sens]
		}
	}
	return(order)
}
