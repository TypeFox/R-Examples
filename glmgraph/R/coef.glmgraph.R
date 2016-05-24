coef.glmgraph <- function(object,lambda1,lambda2,...) {
  if (missing(lambda1) && missing(lambda2)) {
  	  	return(object$betas)
  }
  else if( missing(lambda1) && !missing(lambda2)  ){
  		ind2 <- match(lambda2,object$lambda2)
  		if (any(is.na(ind2))) 
  			stop(paste( "lambda2 should be specified as the following list",object$lambda2))
  		betas <- list()
  		for(i in 1:length(ind2)) betas[[i]] <- object$betas[[ind2[i]]]
  		names(betas) <- lambda2
  		return(betas)
  }
  else if (!missing(lambda1) && missing(lambda2)) { #average lambda1 for each lambda2
   if( any(lambda1 > max(object$lambda1)) || any(lambda1 < min(object$lambda1)) )
  	  	 stop("Specified lambda1 should be within the range of lambda1 used to train the model") 
  	  betas <- list()
  	  for(i in 1:length(object$lambda2)){
  	  	 ind1 <- match(lambda1,object$lambda1s[[i]])
  	  	 if (any(is.na(ind1))){
  	  	 	if( any(lambda1 > max(object$lambda1s[[i]])) || any(lambda1 < min(object$lambda1s[[i]])) ){
  	  	 		warning(paste("Specified lambda1:", lambda1, "should be within the range of lambdal for current lambda2:",lambda2)) 
				next;
  	  	 	}	
  	  	 	ind <- approx(object$lambda1s[[i]],seq(object$lambda1s[[i]]),lambda1)$y
    		l <- floor(ind)
    		r <- ceiling(ind)
    		w <- ind %% 1
    		betas[[i]] <- (1-w)*object$betas[[i]][,l,drop=FALSE] + w*object$betas[[i]][,r,drop=FALSE]
  	  	 }else{
  	  	 	betas[[i]] <- object$betas[[i]][,ind1,drop=FALSE]
  	  	 }
  	  }
  	  idx <- !sapply(betas, is.null)
  	  betas <- betas[idx]
  	  names(betas) <- object$lambda2[idx]
  	  return(betas)  	  
  }
  else if(!missing(lambda1) && !missing(lambda2)){
  		ind2 <- match(lambda2,object$lambda2)
  		if (any(is.na(ind2))) 
  			stop(paste( "lambda2 should be specified as the following list",object$lambda2))
  	    if( any(lambda1 > max(object$lambda1)) || any(lambda1 < min(object$lambda1)) )
  	  	 	stop("Specified lambda1 should be within the range of lambda1 used to train the model") 
  	  	betas <- list()
  	  	for(i in 1:length(lambda2)){
  	  		ind1 <- match(lambda1,object$lambda1s[[ind2[i]]])
  	  	 	if (any(is.na(ind1))){
  	  	 		if( any(lambda1 > max(object$lambda1s[[ind2[i]]])) || any(lambda1 < min(object$lambda1s[[ind2[i]]])) ){
  	  	 			 warning(paste("Specified lambda1:", lambda1, "should be within the range of lambdal for current lambda2:",lambda2[i])) 
					 next;
  	  	 		}
  	  			ind <- approx(object$lambda1s[[ind2[i]]],seq(object$lambda1s[[ind2[i]]]),lambda1)$y
    			l <- floor(ind)
    			r <- ceiling(ind)
    			w <- ind %% 1
    			betas[[i]] <- (1-w)*object$betas[[ind2[i]]][,l,drop=FALSE] + w*object$betas[[ind2[i]]][,r,drop=FALSE]
  	  		}else{
  	  			betas[[i]] <- object$betas[[i]][,ind1,drop=FALSE]
  	  		}
  	  	}  
  	  	idx <- !sapply(betas, is.null)
  	  	betas <- betas[idx] 
  	  	names(betas) <- lambda2[idx]
  	  	return(betas)  
  }
}







