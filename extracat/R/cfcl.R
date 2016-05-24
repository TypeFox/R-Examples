cfcl <- function( x, y = NULL, ll ){
	Freq <- NULL
	if(is.null(y)){
		
		
		x <- as.data.frame(x)
		
		stopifnot(dim(x)[2]>1)
		if("Freq" %in% names(x)){
			Freq <- x$Freq
		}
		
		y <- x[,1]
		x <- x[,2]
		
	}
	
	x <- as.factor(x)
	y <- as.factor(y)
	
	for(i in seq_along(ll[[1]])){
		
		levels(y)[match(ll[[1]][[i]],levels(y)) ] <- paste("RG",i,sep="")
		levels(x)[match(ll[[2]][[i]],levels(x)) ] <- paste("CG",i,sep="")
		
	}
	x <- factor(x, levels=sort(levels(x)))
	y <- factor(y, levels=sort(levels(y)))
	if(!is.null(Freq)){
		return(data.frame(y,x,Freq))
	}
	return(data.frame(y,x))
	
}


rccl <- function(cc,i=1){
	ords <- attr(cc,"orders")[[i]]
	N <- sum(sapply(ords,length))
	ret<-sapply(ords,function(z){
		1:N %in% z
		})
	ret <- apply(ret,1,which)
	return(ret)
}


combcl = function(x){
	
	x <- as.data.frame(sapply(x, as.integer))
	r <- rep(0,nrow(x))
	
	ss <- apply(x,1,function(z) diff(z)[1])
	diags <- which(ss==0)
	r[diags] <- x[diags,1]
	r <- factor(r, levels = c(1:max(r),0))
	
	return(r)
	
}
