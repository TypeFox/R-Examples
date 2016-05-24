clust_diss <-
function(A,B) {
	A <- as.matrix(A)
	B <- as.matrix(B)
	AUB <- cbind(A,B)
	repA <- clusterscore(A)
	repB <- clusterscore(B)
	repAUB <- clusterscore(AUB)
	valproA <- repA$sv^2
	valproB <- repB$sv^2
	valproAUB <- repAUB$sv^2
	crit <- valproA+valproB-valproAUB
	return(crit)
}

