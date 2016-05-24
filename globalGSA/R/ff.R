ff <-
function(data, model="add") {
	Geno <- data[,2:ncol(data)]
    	if(model=="add") {
		Geno <- data.frame(sapply(Geno, Codadd))		
	}
	if(model=="codom") {
		Geno <- data.frame(sapply(Geno, Codcodom))
	}
	if(model=="dom") {
    		Geno <- data.frame(sapply(Geno, Coddom))
	}
	if(model=="rec") {
    		Geno <- data.frame(sapply(Geno, Codrec))
	}
	return(data.frame(data[,1],Geno))
}
