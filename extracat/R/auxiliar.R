# imat <- function(s)
		# {
			# s <- as.factor(s)
			# n <- length(s)
			# s <- as.factor(s)
			# x <- matrix(0, n, length(levels( s )) )
			# x[(1:n) + n*(unclass(s)-1)] <- 1
			# dimnames(x) <- list(names(s), levels(s))
			# return(x)
		# }
		
		
		
Burt = function(x){
	if(inherits(x,"table")){
        return(Burt.table(x))
    }
	if(!"Freq" %in% names(x)){
	   x <- subtable(x, 1:ncol(x))
	}
	   fi <- which(names(x) == "Freq")
	   nd <- ncol(x)-1
	Z <- do.call(cbind, sapply(x[, 1:nd], imat, simplify = FALSE))
	rownames(Z)<- NULL
	Z <- as.data.frame(Z)
	Z <- t(Z * x[, fi]) %*% as.matrix(Z)
	if(length(unique(rownames(Z))) < length(rownames(Z)) ){
		nlvl <- sapply(x[,-fi],function(z) nlevels(as.factor(z)))
		rn <- paste( rep(names(x)[-fi],nlvl), rownames(Z), sep=":")	
		rownames(Z) <- colnames(Z) <- rn
	}
return(Z)
}


 # idat = function(x, allcat = FALSE, freqvar = NULL){
	# if("Freq" %in% names(x)) freqvar <- "Freq"
	
	# if(!is.null(freqvar)){
		# fi <- which(names(x) ==freqvar)
		# s <- x[,fi]
		# x <- x[,-fi]	
	# }else{
		# s <- NULL	
	# }
	# if(allcat){
		# ret<-lapply(x,imat)
	# }else{
		# ret<-lapply(x,function(z){
			# y<-imat(z)
			# y <- y[,-ncol(y),drop=FALSE]
		# })
	# }
	# nlvl <- sapply(ret,ncol)
	
	# ret <- as.data.frame(do.call(cbind,ret))
	
	# names(ret) <- paste(  rep(names(x),nlvl), names(ret),sep = ":")
	# attr(ret,"var.ids") <- rep(1:ncol(x),nlvl)
	
	# if(!is.null(s)){
		# ret[freqvar] <- s	
	# }
	# return(ret)
# }
Burt.table <- function(x){
	stopifnot(inherits(x,"table"))
	nd <- length(dim(x))
B2 <- NULL
for(i in 1:nd){
	B1 <- NULL
	for(j in 1:nd){
		B0<-apply(x,unique(c(i,j)),sum)
        if(i == j){
            B0 <- diag(B0)
        }
        colnames(B0) <- dimnames(x)[[j]]
        rownames(B0) <- dimnames(x)[[i]]
		 B1 <- cbind(B1,B0)
		}
		B2 <- rbind(B2,B1)
} 

if(length(unique(rownames(B2))) < length(rownames(B2)) ){
		
		rn <- paste( rep( names(attr(x,"dimnames")) , dim(x) ), rownames(B2), sep=":")	
		rownames(B2) <- colnames(B2) <- rn
	}

return(B2)
}


isSymMat <- function(x,tol = 1e-12){
	if(length(dim(x)) != 2) return(FALSE)
	
	if(diff(dim(x)) != 0) return(FALSE)
	
	err <- sum(abs(x-t(x)))
	if(err<tol){
		return(TRUE)
	}
	return(FALSE)
}
