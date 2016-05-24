
###################################################
# pooling function for nested multiple imputation objects
pool.mids.nmi <- function( object , method="largesample" ){

		call <- match.call()
		Nimp <- object$Nimp
		NB <- Nimp["between"]
		NW <- Nimp["within"]
		anal <- object$analyses

		pool_results <- as.list( 1:NB)
		for (bb in 1:NB){
		#    bb <- 1
			anal.bb <- anal[[bb]]
			anal2 <- list( "analyses" = anal.bb )
			class(anal2) <- "mira"
			pool_results[[bb]] <- mice::pool( anal2 , method=method)
						}

		pool1 <- pool_results[[1]]$u[1,,]
		
		
		NV <- nrow(pool1)
		cn <- colnames(pool)
		v1 <- vcov( object$analyses[[1]][[1]] )
		if ( ! is.null( rownames(v1) ) ){
			cn <- rownames(v1)
								}
		
		# collect parameter estimates
		qhat <- array( NA , dim=c( NB , NW , NV ) )
		dimnames(qhat)[[3]] <- cn
		dimnames(qhat)[[1]] <- paste0("Between_Imp" , 1:NB ) 
		dimnames(qhat)[[2]] <- paste0("Within_Imp" , 1:NW )
		for (bb in 1:NB){
			qhat[bb,,] <- pool_results[[bb]]$qhat
					}

		# collect estimated variance matrices
		u <- array( NA , dim=c( NB , NW , NV , NV) )
		dimnames(u)[[4]] <- dimnames(u)[[3]] <- cn
		dimnames(u)[[1]] <- paste0("Between_Imp" , 1:NB ) 
		dimnames(u)[[2]] <- paste0("Within_Imp" , 1:NW )

		for (bb in 1:NB){
			u[bb,,,] <- pool_results[[bb]]$u
					}
		
		#*************************
		# summary statistics
		fit <- pool.nmi.scalar.helper( qhat , u , NV , NB , NW )
		fit$Nimp <- Nimp
		class(fit) <- "mipo.nmi"
		return(fit)
			}
###########################################################################			
# coef and vcov method
coef.mipo.nmi <- function( object , ... ){
   object$qbar
		}
vcov.mipo.nmi <- function( object , ... ){
   object$Tm
		}		