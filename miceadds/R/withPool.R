

############################################################
# extraction method for multiply imputed datasets
withPool_MI <- function(x,...){
	res <- x
	if (inherits(res,"mira")){ 
			res <- res$analyses
				}
	NB <- length(res)
	res1 <- res[[1]]	
	if (NB>1){
		for (bb in 2:NB){	
			res1 <- res1 + res[[bb]]				
				}
			}
	res1 <- res1 / NB
	return(res1)			
}
#################################################################



############################################################
# extraction method for multiply imputed datasets
withPool_NMI <- function(x,...){
	res <- x
	if (inherits(res,"mira.nmi")){ 
			res <- res$analyses
				}
	NB <- length(res)
	NW <- length(res[[1]])
	res1 <- 0	
	for (bb in 1:NB){
		for (ww in 1:NW){
			res1 <- res1 + res[[bb]][[ww]]				
				}
			}
	res1 <- res1 / (NB*NW)
	return(res1)			
}
#################################################################