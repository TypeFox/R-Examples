#classcrit = function(M){
#	m = ncol(M)
#	n = nrow(M)
#	
#	M2 = M[1:(n-1),m:2]
#	
#	R = apply(M2,1,csm)
#	R2 = apply(R,1,csm)[,(m-1):1]
#	
#	crit = sum(as.numeric(M[2:n,1:(m-1)])*R2)
#	return(crit)
#}
classcrit <- function(x, concordant = TRUE){

	if(!concordant){
		if(length(dim(x)) > 2){
			cat("discordant criterion currently only available for two-way data. Using concordant pairs instead.")	
		}else{
			x <- x[,ncol(x):1]
		}	
	}	
	storage.mode(x) <- "integer"
	return(.Call("classcrit",x,as.integer(dim(x)),as.integer(0)))
}

csm = function(x){
	cumsum(as.numeric(x))
}	
	
hammcrit = function(x){
	stopifnot(length(dim(x))==2)
	
	m = ncol(x)
	n = nrow(x)
	
	x = t(apply(x,1,rev))
	
	RC = apply(x,1,csm)
	RR = apply(x,2,csm)
	

	
	RC2 = apply(RC,1,csm)-t(RC)
	RR2 = t(apply(RR,1,csm))-RR
	
	RC3 = apply(RC2,2,csm)
	RR3 = apply(RR2,1,csm)
	
	crit = sum(as.numeric(x*(RC3+t(RR3))))
	return(crit)
}






#indep.class.crit = function(M){
#	ttx = rowSums(M)
#	tty = colSums(M)
#	
#	piv = outer(tty, tty)/sum(tty)
#                  piv = sum(piv[upper.tri(piv)])
#    pjv = outer(ttx, ttx)/sum(ttx)
#                  pjv = sum(pjv[upper.tri(pjv)])
#                  worst = piv * pjv
#return(worst)
#}

iccrit = function(x){
	crt <- 1
	n <- sum(x)
	if(n == 0) return(0)
	if(is.null(dim(x))){
	  ss <- x/n
	   piv <- outer(ss,ss)
	   crt <- crt*sum(piv[upper.tri(piv)])
	}else{
	nd <- length(dim(x))
	
	for(i in seq_along(dim(x))){
		ss <- apply(x,i,sum)/n
		piv <- outer(ss,ss)
		crt <- crt*sum(piv[upper.tri(piv)])
	}
	}
	crt <- crt*n^2
	
return(crt)
}	
	
#indepcrit2 = function(M, vs = 1){
#	tty = rowSums(M)
#	ttx = colSums(M)
#	n <- nrow(M)
#	m <- ncol(M)
#		yo <- c(seq(1,n,2),rev(seq(2,n,2)))
#		xo <- c(seq(1,m,2),rev(seq(2,m,2)))
#	if(vs == 1){
#		M <- (M[order(tty), order(ttx)])[yo,xo]
#	}else{
#		M <- (M[order(tty, decreasing = TRUE), order(ttx, decreasing = TRUE)])[yo,xo]
#	}
#	
#	return(fastRcrit2(M))
#}

ihcrit = function(x){
	tty = rowSums(x)
	ttx = colSums(x)
	n <- nrow(x)
	m <- ncol(x)
	if((N<-sum(x)) == 0) return(0)
		yo <- c(seq(1,n,2),rev(seq(2,n,2)))
		xo <- c(seq(1,m,2),rev(seq(2,m,2)))
	x <- (tty %*% t(ttx))/N
	x <- (x[order(tty), order(ttx)])[yo,xo]
	return(hammcrit(x))
}	

WBCI = function(x){
	ih <- ihcrit(x)
	if(abs(ih) < 1e-12){return(1)}
	hammcrit(x)/ih
}


sccrit = function(x, concordant = FALSE){
	if(!concordant){
		if(length(dim(x)) > 2){
			cat("discordant criterion currently only available for two-way data. Using concordant pairs instead.")	
		}else{
			x <- x[,ncol(x):1]
		}	
	}	
	storage.mode(x) <- "integer"
	ret <- .Call("classcrit",x, as.integer(dim(x)),as.integer(0))
	return(ret/iccrit(x))
}



itab <- function(x){
   	dims <- dim(x)
    nd <- length(dims)
    N <- sum(x)
    
    p1 <- apply(x,1,sum)
    p2 <- apply(x,2,sum)
    IM <- p1 %o% p2 /N
    if(nd > 2){
        for(j in 3:nd){
            p3 <- apply(x,j,sum)
            IM <- IM %o% p3/N
        }
    }
    dim(IM) <- dim(x)
    return(IM)
}
indep.table <- ind.table <- itab






neg3t <- function(x){
	
	n <- dim(x)[1]
	m <- dim(x)[2]
	k <- dim(x)[3]
	
#c1 <- classcrit(x)	
	c2 <- classcrit(x[n:1,,]) + classcrit(x[,m:1,]) + classcrit(x[,,k:1]) 
	
	return(c2/3/iccrit(x))	
	
}



neg3j <- function(x,r = 1){
	
	n <- dim(x)[1]
	m <- dim(x)[2]
	k <- dim(x)[3]
	
#c1 <- classcrit(x)	
	c2 <- classcrit(x[n:1,,]) + classcrit(x[,m:1,]) + classcrit(x[,,k:1]) 
	ix <- iccrit(apply(x,r,sum))
	x2 <- apply(x,c(1:3)[-r],sum)
	
	return(c2/ix/( classcrit(x2)+2*classcrit(x2,FALSE) )*sum(x)^2)	
	
}


neg3jb <- function(x,r = 1){
	
	n <- dim(x)[1]
	m <- dim(x)[2]
	k <- dim(x)[3]
	
	x2 <- apply(x,c(1:3)[-r],sum)
	v1<-3*iccrit(x2)/( classcrit(x2)+2*classcrit(x2,FALSE) )*neg3t(x)
	
	return(v1)
	
}


neg3c <- function(x,r = 1){
	
	n <- dim(x)[1]
	m <- dim(x)[2]
	k <- dim(x)[3]
	
	xz <- apply(x,c(c(1:3)[-r][1],r),sum)
	yz <- apply(x,c(c(1:3)[-r][2],r),sum)
	
	bxz <- classcrit(xz)
	byz <- classcrit(yz)
	
	nxz <- bxz + classcrit(xz, FALSE)
	nyz <- byz + classcrit(yz, FALSE)
	
	iz <- iccrit(apply(x,r,sum))
	
	v1 <- (nxz*nyz)/iz - (bxz*byz)/iz
	
	v1<- neg3t(x)*3*iccrit(x)/v1
	
	return(v1)
	
}

allpairs <- function(x){
	
	
	N <- sum(x)
	nd <- length(dim(x))
	nl <- 
	
	nl <- c(N*(N-1),0)
	
	if(nd < 2){
		return(0)
	}
	
	for(i in 1:nd){
		tt <- apply(x,i,sum)
		nl[2] <- nl[2] + sum(tt*(tt-1))
	}
		
	
		nl <- c(nl, sapply(2:nd, function(z){
			   comb <- combn(1:nd,z)
			   sum(apply(comb,2,function(y){
						 tt <- apply(x,y,sum)
						sum( tt*(tt-1) )
				}))
		}))
	
	nls <- sum( nl * (-1)^(2:(nd+2)))/2
	
	return(nls)
	
}
#BCI <- function(x,...){
#    UseMethod("BCI")
#}
#BCI.data.frame <- function(){
#
#}

BCI <- function(x){
	nd <- length(dim(x))
	ic <- iccrit(x)
	if(abs(ic) < 1e-12){return(1)}
	
	ret <- (allpairs(x)-classcrit(x))/ic/(2^(nd-1)-1)
	return(ret)
}



BCC <- function(x){
	nd <- length(dim(x))
	ret <- (allpairs(x)-classcrit(x))
	return(ret)
}




cohen <- function(x){
	stopifnot(length(dim(x))==2)
	N <- sum(x)
	sq<-lapply(dim(x),function(z){
		(0:(z-1))/(z-1)
	})
	D1<-sapply(sq[[1]],function(z){
		abs(z-sq[[2]])
	})
	D2<-sapply(sq[[2]],function(z){
		abs(z-sq[[1]])
	})
	
	D<-1-(t(D1)+D2)/2
	
	dt <- sum(D*x)/N 
	ch <- sum(itab(x)*D)/N
	kappa <- (dt-ch)/(1-ch)
	return(kappa)
	
}


ME <- function(x){
	nd <- length(dim(x))
	Z <- lapply(dim(x),function(z){
		ret <- rbind(cbind(0,diag(1,z-1,z-1)),0)
		diag(ret) <- 1
		ret
	})
	ret <- 0
	for(i in 1:nd){
		M <- apply(x,i,as.vector)
		for(j in 1:(ncol(M)-1)){
			ret <- ret + sum(M[,j]*M[,j+1])
		}
		
	}
	return(ret)
	
}


optME<-
function (x, dims = NULL,nstart = 1, solver = "nearest_insertion", return.table = TRUE, adjust.dist = FALSE) 
{
    nd <- length(dim(x))
    
    ind <- c(1:nd)
    
   if(is.null(dims)){
   	dims <- ind
   }
  
    D <- lapply(ind, function(d) {
    	
    		if(! d %in% dims ){
    			return(1:dim(x)[d])
    		}
    	
        M <- apply(x, d, as.vector)
        d <- t(M) %*% M
        d <- max(d) - d
        
        if(adjust.dist){
    		edist <- dist(t(M))
    		d <- d + as.matrix(edist)/10
    	}
    	
        ts <- TSP(d)
        ts <- insert_dummy(ts, n = 1, label = "cut")
        return(ts)
    })
  

# require(TSP)
    ords <- lapply(D, function(d) {
    	if(!inherits(d, "TSP")){
    		return(d)
    	}
    	st <- system.time({
        ns <- min(nstart, attr(d, "Size"))
        sp <- sample(1:attr(d, "Size"), ns, replace = FALSE)
        tl1 <- Inf
        for (i in sp) {
            sts <- solve_TSP(d, method = solver, control = list(start = as.integer(i)))
            tl0 <- attr(sts, "tour_length")
            if (tl0 < tl1) {
                ord <- cut_tour(sts, "cut")
                tl1 <- tl0
            }
        }
        
        })
        return(ord)
    })
    if(return.table){
    	data <- as.data.frame(as.table(x))
    	for (i in 1:nd) {
       	 data[, i] <- factor(data[, i], levels = levels(as.factor(data[, 
           	 i]))[ords[[i]]])
    	}
    	tt <- xtabs(Freq ~ ., data = data)
    	attr(tt,"orders") <- ords
    	return(tt)
     }
    return(ords)
}



