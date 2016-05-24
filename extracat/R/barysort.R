barysort = function(x, vs = 1){
	n <- nrow(x)
	m <- ncol(x)
	
	

	rs <- rowSums(x)
	cs <- colSums(x)
	ro <- 1:n
	co <- 1:m
	
	optimal <- FALSE
	crt <- classcrit(x)
	while(!optimal){
		
		b0 <- ( x %*% 1:m )/rs
		ord0 <- order(b0)
		x <- x[ord0,]
		rs <- rs[ord0]
		ro <- ro[ord0]
		crt0 <- classcrit(x)
		
		b1 <- ( t(1:n) %*% x )/cs
		ord1 <- order(b1)
		x <- x[,ord1]
		cs <- cs[ord1]
		co <- co[ord1]
		
		crt1 <- classcrit(x)
		
		if(crt1 <=crt){
			#column change not good
			# but accepted if row change was ok
			if(crt0 <=crt){
				# row changes are also bad, break here
				optimal <- TRUE
				x <- x[order(ord0),order(ord1)]
				ro <- ro[order(ord0)]
				co <- co[order(ord1)]	
			}else{
				crt <- crt0
				
			}
		}else{
			crt <- crt1
		}
		
		
		
		
	}
	attr(x,"orders") <- list(ro,co)
	return(x)
}


barysortXD <- function(x){

	nd <- length(dim(x))
	mars <- lapply(dim(x), function(z) 2*(1:z)/z/(z+1))
	optimal <- FALSE
	dx <- as.data.frame(as.table(x))
	Ix <- list()
	for(i in 1:nd){
		tmp <- mars
		#tmp[[i]] <- rep(1/dim(x)[i],dim(x)[i])
		Ix[[i]] <- getIx(tmp)
	}
	print(sapply(Ix,dim))
	w <- list()
	for(i in 1:nd){
		w[[i]] <- apply(x,i,sum)
	}
	ind <- lapply(dim(x), function(z) 1:z)
	
	while(!optimal){
		optimal <- TRUE
		for( i in 1:nd ){
			M <- x * Ix[[i]]
			b <-apply(M,i,sum)/w[[i]]
			ord <- order(b)
			if( !all( ord == 1:length(ord))){
				# accept new order
				dx[,i] <- factor(dx[,i], levels = levels(dx[,i])[ord])
				w[[i]] <- w[[i]][ord]
				ind[[i]] <- ind[[i]][ord]
				x <- xtabs(Freq~., data = dx)
				optimal <- FALSE

			}
			
		print(BCI(x))
			
		}
		
	
	
	}
	attr(x,"ord") <- ind
	return(x)

}

getIx <- function(x){

	nd <- length(x)
	if(nd == 1) return(nd)

	P <- outer(x[[1]],x[[2]])
	if(nd == 2) return(P)
	
	for(i in 3:nd){
		P <- outer(P,x[[i]])
	}
	return(P)
	
}


#presortMat <- function(x, sym = FALSE, v = 1){
#	stopifnot(length(dim(x))==2)
#	n <- nrow(x)
#	m <- ncol(x)
#	if(v == 1){
#		mp1 <- seq((n-1),1-n,-2)
#		mp2 <- seq((m-1),1-m,-2)
#		
#		opt <- FALSE
#		ordx0 <- 1:dim(x)[1]
#		ordy0 <- 1:dim(x)[2]
#		ordx <- ordx0
#		ordy <- ordy0
#		
## row profiles
#		rx <- x/rowSums(x)
#		
##column profiles
#		cx <- x/ rep(colSums(x),each = nrow(x))
#		
#		while(!opt){
#			opt <- TRUE
#			
##column profiles
#			val <- mp1 %*% cx
#			ordx1 <- order(val, decreasing = TRUE)
#						
#			if(!all(ordx1==ordx0)){
##change
#				ordx0 <- ordx1
#				opt <- FALSE
#				cx <- cx[,ordx1]
#				rx <- rx[,ordx1]
#				ordx <- ordx[ordx1]
#			}	
#			
##row profiles anew		
#			val <- rx %*% mp2
#			ordy1 <- order(val,decreasing =TRUE)
#			if(!all(ordy1==ordy0)){
##change
#				ordy0 <- ordy1
#				opt <- FALSE
#				cx <- cx[ordy1,]
#				rx <- rx[ordy1,]
#				ordy <- ordy[ordy1]
#			}
#		}
#		return(x[ordx,ordy])
#		
#	}
#	
#	if(v == 2){
#		
#		opt <- FALSE
#		ordx0 <- 1:dim(x)[1]
#		ordy0 <- 1:dim(x)[2]
#		ordx <- ordx0
#		ordy <- ordy0
#		
#		
## row profiles
#		rx <- x/rowSums(x)
#		
##column profiles
#		cx <- x/ rep(colSums(x),each = nrow(x))
#		
#	
#	while(!opt){
#			opt <- TRUE
#			
##column profiles
#				mp1 <- apply(rx,2,mean)
#			val <- mp1 %*% cx
#			ordx1 <- order(val, decreasing = TRUE)
#						
#			if(!all(ordx1==ordx0)){
##change
#				ordx0 <- ordx1
#				opt <- FALSE
#				cx <- cx[,ordx1]
#				rx <- rx[,ordx1]
#				ordx <- ordx[ordx1]
#			}	
#			
##row profiles anew		
#				mp2 <- apply(cx,1,mean)
#			val <- rx %*% mp2
#			ordy1 <- order(val,decreasing =TRUE)
#			if(!all(ordy1==ordy0)){
##change
#				ordy0 <- ordy1
#				opt <- FALSE
#				cx <- cx[ordy1,]
#				rx <- rx[ordy1,]
#				ordy <- ordy[ordy1]
#			}
#		}
#		return(x[ordx,ordy])
#		
#	}
#	
#}
#
#
