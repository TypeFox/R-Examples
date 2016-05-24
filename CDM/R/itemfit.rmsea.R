
###########################################################
# RMSEA Item fit
itemfit.rmsea <- function( n.ik , pi.k , probs , itemnames=NULL){
	# probs ... [ classes , items , categories ]
	# n.ik ... [ classes , items , categories , groups ]	
	# N.ik ... [ classes , items , categories]		

	# RMSEA for all groups
	itemfit.rmsea <- .rmsea.aux( n.ik , pi.k , probs )
	if ( ! is.null(itemnames) ){
		names(itemfit.rmsea) <- itemnames }
	# groupwise RMSEA
	G <- dim(n.ik)[4]
	I <- dim(n.ik)[2]
	rmsea.groups <- matrix( NA , I , G )
	if ( ! is.null(itemnames) ){
		rownames(rmsea.groups) <- itemnames }
	for (gg in 1:G){
		rmsea.groups[,gg] <- .rmsea.aux( n.ik[,,,gg,drop=FALSE] , 
					pi.k , probs )
						}	
	res <- list( "rmsea" = itemfit.rmsea , 
				 "rmsea.groups"=rmsea.groups )
	return(res)
	}
##########################################
	
##########################################
# auxiliary function
.rmsea.aux <- function( n.ik , pi.k , probs , eps=10^(-30) ){
	# probs ... [ classes , items , categories ]
	# n.ik ... [ classes , items , categories , groups ]	
	# N.ik ... [ classes , items , categories]	
	N.ik <- n.ik[,,,1]
	G <- dim(n.ik)[4]
	pitot <- pi.k[,1]
	eps <- 10^(-10)
	if (G>1){ 
		for (gg in 2:G ){
			N.ik <- N.ik + n.ik[,,,gg]
			pitot <- pitot + pi.k[,gg]
				}
			}
			
	# calculate summed counts
	N.ik_tot <- array( 0 , dim=dim(N.ik) )
	N.ik_tot[,,1] <- N.ik[,,1,drop=FALSE]
	K <- dim(N.ik)[3]			
	for (kk in 2:K){
		N.ik_tot[,,1] <- N.ik_tot[,,1,drop=FALSE] + N.ik[,,kk,drop=FALSE] 
					}

	for (kk in 2:K){	N.ik_tot[,,kk] <- N.ik_tot[,,1] }
	
	# calculate itemwise statistics
	p.ik_observed <- N.ik / ( N.ik_tot + eps )
	p.ik_observed[ is.na( p.ik_observed ) ] <- 0
	
	# define class weights 
	pi.k_tot <- array( 0 , dim=dim(p.ik_observed) )
	for (kk in 1:K){
#		pi.k_tot[,,kk] <- matrix( pitot , nrow= dim(pi.k_tot)[1] , ncol=dim(pi.k_tot)[2] , byrow=T )
		pi.k_tot[,,kk] <- matrix( pitot , nrow= dim(pi.k_tot)[1] , ncol=dim(pi.k_tot)[2] , byrow=FALSE )
				}
	# calculate statistics
	dist.item <- pi.k_tot * ( p.ik_observed - probs )^2			
# Revalpr("round(pi.k_tot[,6,2],4)")	
# Revalpr("round(p.ik_observed[,6,2],5)")
# Revalpr("round(probs[,6,2],5)")
# h1 <- cbind( pi.k_tot[,6,2] , p.ik_observed[,6,2] , probs[,6,2] )
# write.csv2( h1 , "mod1.csv" )
	h1 <- dist.item[,,1]
	for (kk in 2:K){ h1 <- h1 + dist.item[,,kk] }
	itemfit.rmsea <- sqrt( colSums( h1 ) )
	return(itemfit.rmsea)
		}
