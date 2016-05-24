subtree = function(tree, k = NULL, h = NULL){
	
	data <- tree$order 
	itr <- untree(tree$merge, tree$order, ind = TRUE)
	
	if( is.null(k) & is.null(h) ){
		stop("Neither k nor h specified")	
	}
	if( !is.null(k) & !is.null(h) ){
		stop("Unclear specification: found both k and h.")	
	}
	if(!is.null(h)){
		k <- sum(tree$height > h)+1	
	}
	
	nit <- length(itr)
	
	itr <- itr[ (nit-k+2):nit ]
	mrg <- tree$merge[(nit-k+2):nit,,drop=FALSE]
	
	#compute level order
	
	tmrg <- t(mrg)
	mins <- sapply(itr,function(x) sapply(x, min))
	#mins[tmrg < (nit-k+2)] <- rank(mins[tmrg < (nit-k+2)])
	#tmrg[tmrg < (nit-k+2)] <- -c(1:k)
	#ord <- rank(mins[tmrg < (nit-k+2)])
	tmrg[tmrg < (nit-k+2)] <- -rank(mins[tmrg < (nit-k+2)])
	tmrg[tmrg >= (nit-k+2)] <- tmrg[tmrg >= (nit-k+2)] - nit + k - 1
	mrg <- t(tmrg)
	rm(tmrg)
	
	# TODO: better solution than a loop?
	#s <- 1
	for( i in 1:(k-1) ){
		i1 <- mrg[i,1]
		if( i1 < 0 ){
			data[ itr[[i]][[1]] ] <- -i1 #-s
			#s <- s+1
		}
		i2 <- mrg[i,2]
		if( i2 < 0 ){
			data[ itr[[i]][[2]] ] <- -i2 #-s
			#s <- s+1
		}
	}
	# CHECK
	labs <- as.vector(t(mrg))
	labs <- labs[labs < 0]
	labs <- match(1:k, -labs)
	# CHECK
	
	ret <- new("list")
	ret$merge <- mrg
	ret$order <- 1:k
	ret$labels <- labs
	#ret$data <- factor(data[ order(tree$order) ], levels <- ord)
	ret$data <- factor(data[ order(tree$order) ])
	ret$height <- tree$height[ (nit-k+2):nit ]
	class(ret) <- "hclust"
	return(ret)
}




noderange <- function(mrg,ord, xval = NULL, ind = FALSE){
	stopifnot(ncol(mrg) == 2)
	ng <- nrow(mrg)
	
	ind1 <- lapply(1:ng,function(z) c(ng+2,0))
	ind2 <- ind1
	for(i in 1:ng){
		if(mrg[i,1] < 0){
			ind1[[i]][2] <- max( match(abs(mrg[i,1]),ord), ind1[[i]][2] )
			ind1[[i]][1] <- min( match(abs(mrg[i,1]),ord), ind1[[i]][1] )
		}else{
			ind1[[i]][1] <- min( ind1[[i]][1] , ind1[[ mrg[i,1] ]][1], ind2[[ mrg[i,1] ]][1])
			ind1[[i]][2] <- max( ind1[[i]][2] , ind1[[ mrg[i,1] ]][2], ind2[[ mrg[i,1] ]][2])
		}
		if(mrg[i,2] < 0){
			ind2[[i]][2] <- max( match(abs(mrg[i,2]),ord), ind2[[i]][2] )
			ind2[[i]][1] <- min( match(abs(mrg[i,2]),ord), ind2[[i]][1] )
		}else{
			ind2[[i]][1] <- min( ind2[[i]][1] , ind1[[ mrg[i,2] ]][1], ind2[[ mrg[i,2] ]][1])
			ind2[[i]][2] <- max( ind2[[i]][2] , ind1[[ mrg[i,2] ]][2], ind2[[ mrg[i,2] ]][2])
		}
	}	
	
	ret <- mapply(function(x,y) list( x,y ), x = ind1,y=ind2, SIMPLIFY =FALSE)	
	return(ret)
}

nodexc <- function(mrg,ord, xval = NULL, ind = FALSE){
	stopifnot(ncol(mrg) == 2)
	ng <- nrow(mrg)

	cnt <- as.list(rep(0,ng))
	
	for(i in 1:ng){
		if(mrg[i,1] < 0){
			a <-  (match(abs(mrg[i,1]),ord)-0.5)/(ng+1)
		}else{
			a <- cnt[[mrg[i,1]]][3]	
		}
		if(mrg[i,2] < 0){
			b <-  (match(abs(mrg[i,2]),ord)-0.5)/(ng+1)
		}else{
			b <- cnt[[mrg[i,2]]][3]	
		}
		cnt[[i]] <- c(a,b,(a+b)/2)
	}	

	return(cnt)
}



nodeheight <- function(mrg,ord,ht){
	stopifnot(ncol(mrg) == 2)
	ng <- nrow(mrg)
	#ind1 <- vector(mode="list",length=ng)
	#ind2 <- vector(mode="list",length=ng)
	LH <- as.list(ht)
	RH <- LH
	for(i in 1:ng){
		if(mrg[i,1] > 0){
			LH[[i]] <- LH[[i]] - ht[[mrg[i,1]]]
		}
		if(mrg[i,2] > 0){
			RH[[i]] <- RH[[i]] - ht[[mrg[i,2]]]
		}
	}	

	ret <- mapply(function(x,y) list( x,y ), x = LH,y=RH, SIMPLIFY =FALSE)	
	return(ret)
}
