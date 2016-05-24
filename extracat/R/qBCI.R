

qBCI <- function(x,...){
	UseMethod("qBCI")
}

qBCI.default <- function(x,y,p = NULL, k = 5,iter=20, ...){
	N <- length(x)
	
	if(is.null(p)){
		p <- sqrt(k/N)
		p <- 1/floor(1/p)
	}
	if(!is.factor(x)){
		qx <- quantile(x,seq(0,1,p),na.rm = TRUE)
		qx[1] <- -Inf
		if( length(qxx <- unique(qx)) < length(qx)){
			simpleWarning("non-unique quantiles detected")	
		}
		if(length(qxx) == 2){
			x <- as.factor(x)	
		}else{
			x <- cut(x, unique(qxx) )
		}
	}
	if(!is.factor(y)){
		qy <- quantile(y,seq(0,1,p),na.rm = TRUE)
		qy[1] <- -Inf
		if( length(qyy <- unique(qy)) < length(qy)){
			simpleWarning("non-unique quantiles detected")	
		}
		if(length(qyy) == 2){
			y <- as.factor(y)	
		}else{
			y <- cut(y, unique(qyy) )
		}
	}
	tt <- table(x,y)
	BCI(optile(tt,iter=iter))
}



qBCI.data.matrix <- function(x,p = NULL, k = 5, sort = TRUE, iter=20, ...){
	
	N <- nrow(x)
	nd <- ncol(x)
	
	
	if(is.null(p)){
		p <- (k/N)^(1/nd)
		p <- 1/floor(1/p)
	}
	
	x <- sapply(x, function(z) {
		qz <- quantile(z,seq(0,1,p),na.rm = TRUE)
		qz[1] <- -Inf
		cut(z, qz)
		})
	x <- subtable(x,1:nd)
	x <- xtabs(Freq~.,data=x)
	if(sort){
		return(BCI(optile(x,iter=iter)))
	}else{
		return(BCI(x))
	}
	
}

qBCI.data.frame <- function(x,p = NULL, k = 5, sort = TRUE, iter=20, ...){
	x <- data.matrix(x)
	NextMethod("qBCI",x = x, p = p, k = k, sort = sort, iter = iter)
}


BCImat <- function(x, k = 5, iter = 20, p = NULL){

	nd <- ncol(x)
	ids <- combn(1:nd,2)
	
	values <- apply(ids,2,function(id){
		z <- na.omit(x[,id])
		qBCI(z[,1],z[,2], k = k, p = p, iter = iter)
	})
	M <- matrix(0,nd,nd)
	M[lower.tri(M)] <- values
	M <- M + t(M)
	colnames(M) <- rownames(M) <- names(x)
		return(M)
}


WBCImat <- function(x, iter = 20, freqvar = NULL, fun = "WBCC"){

	if("Freq" %in% names(x)) freqvar <- "Freq"
	
	if(!is.null(freqvar)){
		fi <- which(names(x)==freqvar)
		names(x)[fi] <- "Freq"	
	}else{
		x <- subtable(x,1:ncol(x))	
		fi <- ncol(x)
	}
	
	nd <- ncol(x)
	ids <- combn( (1:nd)[-fi],2)
	
	values <- apply(ids,2,function(z){
		ret<-WBCI(tt<-optile(xtabs(Freq~x[,z[1]]+x[,z[2]],data=x),iter=iter, fun = fun))
		return(ret)
	})
	M <- matrix(0,nd-1,nd-1)
	M[lower.tri(M)] <- values
	M <- M + t(M)
	colnames(M) <- rownames(M) <- names(x)[-fi]
		return(M)
}

#BCImat2 <- function(x, k = 5, iter = 20, p = NULL){
	
#	nd <- ncol(x)
#	ids <- combn(1:nd,2)
	
#	types <- sapply(x, typeof)
#	if(any(types != "factor")){
#		vi <- which(types != "factor")
#		x[,vi] <- sapply(x[,vi], function(z) qbin(z, k = k, p = p))
#	}
	
	
#	values <- apply(ids,2,function(id){
#					tt <- table(x[,id[1]],x[,id[2]])
#					BCI(optile(tt, iter = iter))
#					})
#	M <- matrix(0,nd,nd)
#	M[lower.tri(M)] <- values
#	M <- M + t(M)
#	colnames(M) <- rownames(M) <- names(x)
#	return(M)
#}



qbin <- function(x,p = NULL, k = 5, d = 2){
	N <- length(x)
	if(is.null(p)){
		p <- (k/N)^(1/d)
		p <- 1/floor(1/p)
	}
	if(!is.factor(x)){
		qx <- quantile(x,seq(0,1,p), na.rm = TRUE)
		qx[1] <- -Inf
		if( length(qxx <- unique(qx)) < length(qx)){
			simpleWarning("non-unique quantiles detected")	
		}
		if(length(qxx) == 2){
			x <- as.factor(x)	
		}else{
			x <- cut(x, unique(qxx) )
		}
	}
}


cmat <- function(x, sort = TRUE, crit = BCI, k = 5, iter = 20, p = NULL, jitter = TRUE, freqvar = NULL, diag = NULL, fun = "BCC", foreign = NULL){

	nd <- ncol(x)
	
	if(!is.null(freqvar)){
		fi <- which( names(x)==freqvar )
		stopifnot(length(fi)>0)
		wts <- x[,fi] 
		x <- x[, -fi]
	}else{
		wts <- rep(1,nrow(x))
	}
	
	factors <- sapply(x, is.factor)
	N <- nrow(x)
	if(!all(factors)){
			if(is.null(p)){
				
				p <- sqrt(k/N)
				p <- 1/floor(1/p)
			}
			
			for(i in which(!factors)){
				xx <- x[,i]
				if(jitter){
					md <- diff(sort(na.omit(xx)))
					md <- min(md[md>0])
					
					xx <- xx + runif(N,-md/2,md/2)
				}
				qx <- quantile(xx,seq(0,1,p),na.rm = TRUE)
				qx[1] <- -Inf
				qxx <- unique(qx)
				
				if( length(qxx) < length(qx)){
					simpleWarning(paste("non-unique quantiles detected in variable ",names(x)[i]))	
				}
				if(length(qxx) == 2){
						x[,i] <- as.factor(x[,i])	
				}else{
					x[,i] <- cut(xx, unique(qxx) )
				}
			}
	}
	
	x$Freq <- wts
	inc <- ifelse(is.null(diag),0,1)
		M <- matrix(0,nd,nd)
		for(i in 1:(nd-1)){
			for(j in (i+inc):nd){
				if(sort){
					M[i,j] <- M[j,i] <- crit(optile(x[,c(i,j,ncol(x))], iter = iter, fun = fun, 
						foreign = foreign,return.type = "table")) 
										
					#M[i,j] <- M[j,i] <- crit(optile(xtabs(wts~x[,i]+x[,j]), iter = iter, fun = fun, foreign = foreign)) 
					
					##attr(optile(table(x[,i],x[,j]), iter = iter) ,"scaled.criterion")
				}else{
					M[i,j] <- M[j,i] <-  crit(xtabs(wts~x[,i]+x[,j]))
				}
				
			}
		}	
	
		if(!is.null(diag)) diag(M) <- diag


	
		colnames(M) <- rownames(M) <- names(x)[1:nd]
		return(M)
}



binary.bci <- function(x, sort = TRUE, crit = BCI, freqvar = NULL, diag = NULL){

	nd <- ncol(x)
	ids <- combn(1:nd,2)
	
	if(!is.null(freqvar)){
		fi <- which( names(x)==freqvar )
		stopifnot(length(fi)>0)
		wts <- x[,fi] 
		x <- x[, -fi]
	}else{
		wts <- rep(1,nrow(x))
	}
	inc <- ifelse(is.null(diag),0,1)
	
	tb <- lapply(x,table)
	N <- nrow(x)
	
		M <- matrix(0,nd,nd)
		for(i in 1:(nd-1)){
			for(j in (i+inc):nd){
				if(sort){
					tt <- xtabs(wts~x[,i]+x[,j])
					M[i,j] <- M[j,i] <- min( tt[1,1]*tt[2,2], tt[1,2]*tt[2,1] )/(prod(tb[[i]])*prod(tb[[j]]))*N^2
					 
				}else{
					M[i,j] <- M[j,i] <-  tt[1,2]*tt[2,1] /(prod(tb[[i]])*prod(tb[[j]]))*N^2
				}
				
			}
		}	
	
	if(!is.null(diag)) diag(M) <- diag


	
		colnames(M) <- rownames(M) <- names(x)
		return(M)
}


