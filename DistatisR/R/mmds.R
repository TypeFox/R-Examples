mmds <-
function(DistanceMatrix,masses=NULL){
	# mmds (metric or classical mds)
	# Quick and dirty mds with or without masses
	# here to accomagny distatis and allow chi2 distance analyses
	# WARNING mds assumes that the distance matrix
	# is (generalized) Euclidean
	# and so it keeps only the positive eigenvalues
	D <- DistanceMatrix   # Being lazy!
	nI <- nrow(D)
	if (is.null(masses)) {masses <- rep(1/nI,nI)} # create masses if needed
	m <- masses # lazy again
	LeM = t(kronecker(m,t(rep(1,nI))))  # PM repmat
	Xhi <- diag(1,nI) -  LeM # centering matrix
	#print(Xhi)                       
	S <- -.5*sqrt(LeM)*Xhi%*%D%*%t(Xhi)*sqrt(t(LeM))   # Double Centered SCP matrix
	#print(S)
	#pause()
	eig <- eigen(S,symmetric=TRUE) # Eigen-decomposition of S
	
	# clean to keep positive eigenvalues only
	Cleaner =  which(eig$value > 0)
	U <- eig$vector[,Cleaner]
    L <- eig$values[Cleaner]
     Nom2Dim = paste('dim',1:length(L))
         names(L) <- Nom2Dim
    Nom2Row = rownames(D)
    # Factor scores
	LeF <- kronecker(1/sqrt(m),t(rep(1,length(L))))*t(t(U) *  sqrt(L) ) 
	
	rownames(LeF) <- Nom2Row
	colnames(LeF) <- Nom2Dim
	# Percentage of inertia
	tau <- round(100*(L/sum(L)),digits=2)
	names(tau) <- Nom2Dim
	# Contributions
	Ctr <- kronecker(m,t(rep(1,length(L))))*t(t(LeF^2) /L )
	rownames(Ctr) <- Nom2Row
	colnames(Ctr) <- Nom2Dim
	#	return(list(S=S,fi=LeF,pdq=eig,M=LeM))
	return(list(FactorScores=LeF,eigenvalues=L,Contributions=Ctr,percentage=tau,M=LeM))
}
