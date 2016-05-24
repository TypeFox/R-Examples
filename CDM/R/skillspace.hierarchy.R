###############################################################
# computation of skill space hierarchy
skillspace.hierarchy <- function( B , skill.names  ){
	
	if ( ! is.matrix(B) ){
		Blist <- strsplit(B , split="\\n")[[1]]
		BL <- length(Blist)

		K <- length(skill.names)
		B <- matrix( 0 , nrow=K , ncol=K)
		rownames(B) <- colnames(B) <- skill.names
		for (bb in 1:BL){
			# bb <- 1
			Blist.bb <- gsub( " " , "" , Blist[[bb]] )
			Bl2 <- strsplit( Blist.bb , split=">" )[[1]]
			VV <- length(Bl2)
			for ( vv in 1:(VV-1) ){
				B[ Bl2[vv] , Bl2[vv+1] ] <- 1
								}
						}
				}		
#	if ( is.null(skill.names) ){  
#		if ( is.null( colnames(B) ) ){
#			skill.names <- paste0( "A" , 1:ncol(B) ) } else {
#				skill.names <- colnames(B) 
#							}
#					}	
	K <- length(skill.names)
	# define complete skill space
	dfr <- rbind( rep(0,K) , rep(1,K) )
	skillspace <- expand.grid( as.list(as.data.frame(dfr) ) )
	colnames(skillspace) <- skill.names
	# attribute pattern labels
	n1 <- paste0("A" , skillspace[,1] )
	for (nn in 2:K){
		n1 <- paste0( n1 , skillspace[,nn] )
				}
	rownames(skillspace) <- n1	
	skillspace0 <- skillspace

	# compute reachability
	R <- B
	V1 <- R
	vv <- 0
	while( ( abs( sum(R) ) > 0 ) & ( vv < 200 ) ){
		R <- R %*% B
		V1 <- V1 + R
		vv <- vv+1
				}
	# exclude skill classes
	for (ii in 1:K){
		for (jj in 1:K){
			if ( ( V1[ii,jj] > 0 ) & ( ii != jj)  ){
				# ii <- 1 ; jj <- 2
				ind <- which( ( skillspace[ , ii ] == 0 ) & ( skillspace[,jj] == 1 ) )
				if ( length(ind) > 0 ){ skillspace <- skillspace[ - ind , ] }
										}
							}
			}
	#**** determine skill classes which were removed	
	zeroprob.skillclasses <- which( ! ( rownames(skillspace0)  %in% rownames(skillspace) ) )
	
	#**************************************
	# output        
	res <- list("R" = V1 , "skillspace.reduced" = as.matrix(skillspace) ,
			"skillspace.complete" = as.matrix(skillspace0) , 
			"zeroprob.skillclasses" =zeroprob.skillclasses )
	return(res)
		}
############################################################################		

# full skill space
skillspace.full <- function( skill.names ){
    B <- paste0( skill.names[1] , " > " , skill.names[2] )
    M1 <- skillspace.hierarchy( B = B , skill.names=skill.names )
    return(M1$skillspace.complete)
                }