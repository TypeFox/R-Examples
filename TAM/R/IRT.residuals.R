

#####################################################
# IRT.residuals
IRT.residuals <- function (object, ...) {
    UseMethod("IRT.residuals")
       }
#####################################################


#######################################################################
tam.residuals <- function( object , ... ){ 
    tamobj <- object
	res <- tam.wle( tamobj , progress=FALSE , 
				output.prob=TRUE , ...  )
	probs <- res$probs
	probs[ is.na(probs) ] <- 0
	theta <- res$theta 
	resp <- tamobj$resp
	B <- tamobj$B
	if ( dim(B)[3] > 1 ){
		stop("Residuals are only computed for unidimensional models.\n")
						}
	B <- B[,,1]	
	K <- dim(B)[2]
	I <- dim(B)[1]
	N <- nrow(resp)
	# compute expected value
    X_exp <- matrix( 0 , nrow=N , ncol=I)
	colnames(X_exp) <- colnames(resp)
	X_exp[ is.na(resp) ] <- NA
	X_var <- X_exp
	for (kk in 1:K){
	    B.kk <- B[,kk]
		B.kk <- matrix( B.kk , nrow=N , ncol=I , byrow=TRUE )
		X_exp <- X_exp + B.kk * aperm( probs[ , kk , ] , c(2,1) )
		X_var <- X_var + B.kk^2 * aperm( probs[ , kk , ] , c(2,1) )
					}
	X_var <- X_var - X_exp^2				
	# compute residuals
	residuals1 <- resp - X_exp
	stand_residuals <- ( resp - X_exp ) / sqrt( X_var )
	# output				
	res <- list( residuals=residuals1 , stand_residuals=stand_residuals , 
				probs=probs , X_exp = X_exp , 
				X_var = X_var , theta=theta , probs = probs)				
	return(res)
		}
#######################################################################		

residuals.tam.mml <- IRT.residuals.tam.mml <- tam.residuals
residuals.tam.mml.2pl <- IRT.residuals.tam.mml.2pl <- tam.residuals
residuals.tam.mml.mfr <- IRT.residuals.tam.mml.mfr <- tam.residuals