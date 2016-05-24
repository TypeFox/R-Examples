makeDistancesAndWeights <- function(DATA,method="euclidean",masses=NULL){
	if(method=="chi2"){
		chi2res <- chi2Dist(DATA)
		D <- chi2res$D
		MW <- list(M=chi2res$M)
	}else if(method=="correlation"){
		#why the hell does this not need to be squared?
		D <- (1-cor(t(DATA)))
		MW <- computeMW(D,masses=masses)
	}else{
		D <- as.matrix(dist(DATA,method=method,diag=TRUE,upper=TRUE))
		D <- D^2
		MW <- computeMW(D,masses=masses)		
	}
	return(list(D=D,MW=MW))
}