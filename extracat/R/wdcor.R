
dcorOld = function(x){
	n <- nrow(x)
	m <- ncol(x)
	
	cases <- expand.grid(1:n,1:m)
	DMY <- as.matrix(dist(cases[,1]))
	DMX <- as.matrix(dist(cases[,2]))
	
	pv <- as.vector(x)/sum(x)
	
	S1 <- pv %*% (DMX*DMY) %*% pv
	S2 <- (pv %*% DMX %*% pv) * (pv %*% DMY %*% pv)
		v1 <- DMY %*% pv
		v2 <- DMX %*% pv
	S3 <- t(v1*v2) %*% pv
	S1X <- pv %*% (DMX*DMX) %*% pv
	S2X <- (pv %*% DMX %*% pv)^2
	S3X <- t(v2^2) %*% pv
	
	S1Y <- pv %*% (DMY*DMY) %*% pv
	S2Y <- (pv %*% DMY %*% pv)^2
	S3Y <- t(v1^2) %*% pv
	
	dcor <- (S1+S2-2*S3)/sqrt(  (S1X+S2X-2*S3X) * (S1Y+S2Y-2*S3Y) )
	return(sqrt(dcor))
}
wdcor = function(x,...){
	UseMethod("wdcor")
}

wdcor.default = function(x,y,w = NULL, ep = 1, approx = FALSE, n = 50, na.rm = TRUE, ...){
	if(length(x) > 20000){
	 return(approx.dcor(x,y,n,ep=ep))
	}
	if(is.null(w)){
		w <- rep(1,length(x))
	}
	stopifnot(length(x) == length(w))
	stopifnot(length(x) == length(y))
	stopifnot(all(w >= 0))
	if(na.rm){
		na <- is.na(x) | is.na(y) | is.na(w)
		if(any(na)){
			ind <- which(!na)
			x <- x[ind]
			y <- y[ind]
			w <- w[ind]
		}
		
	}
	
	
	storage.mode(w) <- "double"
	storage.mode(x) <- "double"
	storage.mode(y) <- "double"
	storage.mode(ep) <- "double"

ret <- .Call("dcorR",x,y,w/sum(w),ep)
	
	return(ret)
	
}
wdcor.table = function(x, ep = 1, ...){
	stopifnot( length(dim(x)) == 2) 
	#dx <- as.data.frame(x)
	#dx <- dx[dx$Freq > 0,]
	dx <- subtable(x,1:2)
	dx <- sapply(dx,as.numeric)
	NextMethod("wdcor", x = dx[,1], y = dx[,2], w = dx[,3], ep = ep)
}


approx.dcor2 = function(x,y, n = 50, correct = FALSE, ep = 1){
	stopifnot(length(x) == length(y))
	x <- cut(x,n)
	y <- cut(y,n)
	z <- table(x,y)
	
	ret <- wdcor(z, ep = ep)
	if(correct) ret <- ret*(1 + 2*sqrt(8)/n^2)
	return(ret)
}


approx.dcor = function(x,y, n = 50, ep = 1, bin = "eq"){
	if(length(x) != length(y)){
		stop("Vector lengths differ.")
	}
	
	if(bin %in% c("e","eq","equ","equi")){
		rx <- range(x,na.rm=TRUE)	
		ry <- range(y,na.rm=TRUE)
		xbr <- seq(rx[1],rx[2]+1e-12, (diff(rx)+1e-12)/n)
		ybr <- seq(ry[1],ry[2]+1e-12, (diff(ry)+1e-12)/n)
	}
	if(bin %in% c("q","quant","quantile")){
		xbr <- unique(quantile(x,seq(0,1,1/n)))
		ybr <- unique(quantile(y,seq(0,1,1/n)))
		xbr[length(xbr)]<-xbr[length(xbr)]+1e-12
		ybr[length(xbr)]<-ybr[length(xbr)]+1e-12
	}

	cx <- cut(x,breaks=xbr,include.lowest=TRUE)
	cy <- cut(y,breaks=ybr,include.lowest=TRUE)

	vx <- tapply(x,cx,mean)
	vy <- tapply(y,cy,mean)

	if(any(is.na(vx))){
		ii <- which(is.na(vx))
		vx[ii] <- (xbr[-length(xbr)]+diff(xbr)/2)[ii]
	}
	if(any(is.na(vy))){
		ii <- which(is.na(vy))
		vy[ii] <- (ybr[-length(ybr)]+diff(ybr)/2)[ii]
	}

	z <- subtable(cbind(cx,cy),1:2)

	ret <- wdcor(vx[z[,1]],vy[z[,2]],w = z[,3],ep=ep)
	return(ret)
}

dcorMV <- function(dx,dy,w = NULL){
	stopifnot(inherits(dx,"dist"))
	stopifnot(inherits(dy,"dist"))
	dx <- as.matrix(dx)
	dy <- as.matrix(dy)
	stopifnot(nrow(dx) == nrow(dy))
	
	if(is.null(w)){
		w <- rep(1,length(nrow(dx)))
	}
	storage.mode(w) <- "double"
	storage.mode(dx) <- "double"
	storage.mode(dy) <- "double"
	ret <- .Call("dcorD",dx,dy,w/sum(w))
}


dcorMVtable <- function(x,ind = 1, method = "euclidean"){
	nd <- length(dim(x))
	data <- subtable(x,1:nd,allfactor=TRUE)
	data.matrix <- sapply(data[,-ncol(data)],as.integer)
	dx <- dist(data.matrix[,ind],method = method)
	dy <- dist(data.matrix[,-ind],method = method)
	w <- data$Freq
	ret <- dcorMV(dx,dy,w)
	return(ret)
}


dcorMVdata <- function(x,ind = 1, method = "euclidean", approx = FALSE){
	x <- as.data.frame(x)
	if(approx != FALSE){
		if(approx == TRUE){
			n <- ceiling(100/ncol(x))
			cat("approx = ",n)
		}else{
			n <- approx
		}
		nd <- ncol(x)
		for(i in 1:nd){
			x[,i] <- cut(x[,i],breaks = n)
		}
		
	
		data <- subtable(x,1:nd)
		w <- data$Freq
		data <- sapply(data[,-ncol(data)],as.integer)
	}else{
		data <- subtable(x,1:ncol(x))
		w <- data$Freq
		data <- data[,-ncol(data)]
	}
	dx <- dist(data[,ind],method = method)
	dy <- dist(data[,-ind],method = method)
	
	ret <- dcorMV(dx,dy,w)
	return(ret)
}


opt.wdcor = function(env,xi,yi){
	return(wdcor(as.table(env$mat[yi,xi])))
}

distcor = function(data, dims, perm.cat, ... ){
	nd <- length(dims)
	stopifnot(nd == 2)
	optimal <- FALSE
	globalbest <- -1
	currcrit <- -1
	bestcrit <- -1
	ci <- 1:dims[2]
	ri <- 1:dims[1]
	tri <- ri
	tci <- ci
	e1 <- new.env()
	e1$mat <- xtabs(Freq~.,data=data)
	#print(e1$mat)
	while(!optimal){
		
		#columns
		for(i in 1:dims[2]){
			for(j in 1:dims[2]){
				if(i < j){
					tci[i:(j-1)] <- tci[(i+1):j]
					tci[j] <- ci[i]
					currcrit <- opt.wdcor(e1, tci, ri)
				}
				if(i > j){
					tci[(j+1):i] <- tci[j:(i-1)]
					tci[j] <- ci[i]
					currcrit <- opt.wdcor(e1, tci, ri)
				}
				if(currcrit > bestcrit){
					bestcrit <- currcrit
					ci <- tci	
				}else{
					tci <- ci	
				}	
			}	
		}
		#rows
		for(i in 1:dims[1]){
			for(j in 1:dims[1]){
				if(i < j){
					tri[i:(j-1)] <- tri[(i+1):j]
					tri[j] <- ri[i]
					currcrit <- opt.wdcor(e1, ci, tri)
				}
				if(i > j){
					tri[(j+1):i] <- tri[j:(i-1)]
					tri[j] <- ri[i]
					currcrit <- opt.wdcor(e1, ci, tri)
				}
				if(currcrit > bestcrit){
					bestcrit <- currcrit
					ri <- tri	
				}else{
					tri <- ri	
				}	
			}	
		}
		if(bestcrit <= globalbest){
			optimal	<- TRUE
		}else{
			globalbest <- bestcrit	
		}
	}
    if(kendalls(e1$mat[ri,ci])< 0){
        ri <- rev(ri)
    }
	return(c(ri-1,ci-1,globalbest))
}

dcov2 = function(x){
	n <- nrow(x)
	stopifnot(ncol(x)==2)
	
	MY <- as.matrix(dist(x[,1]))
	MX <- as.matrix(dist(x[,2]))
	
	xc <- apply(MX,2,mean)
	xr <- apply(MX,1,mean)
	xm <- mean(MX)
	
	yc <- apply(MY,2,mean)
	yr <- apply(MY,1,mean)
	ym <- mean(MY)
	
	MX <- (MX - xr) - rep(xc, each=n) + xm
	MY <- (MY - yr) - rep(yc, each=n) + ym
	 dcov <- sum(MX*MY)/n^2
	return(dcov)	
}

#is.wholenumber <- function(x, tol = .Machine$double.eps^0.5){
#	abs(x - round(x)) < tol
#}
#	if(all(is.wholenumber(x))){
#		storage.mode(x) <- "integer"
#	}


# wdcor.data.frame <- function(x, approx = TRUE, n = 50, ...){
	# nmz <- names(x)
	# x <- data.matrix(x)
	# nd <- ncol(x)
	# ids <- combn(1:nd,2)
	# if(approx){
	# values <- apply(ids,2,function(id){
		# z <- na.omit(x[,id])
		# approx.dcor(z[,1],z[,2], n=n)
	# })
	# }else{
		# values <- apply(ids,2,function(id){
		# z <- na.omit(x[,id])
		# wdcor(z[,1],z[,2], n=n)
	# })
	# }
	# M <- matrix(0,nd,nd)
	# M[lower.tri(M)] <- values
	# M <- M + t(M)
    # diag(M) <- 1
	# colnames(M) <- rownames(M) <- nmz
		# return(M)
# }



# wdcor.data.frame <- function(x, w = NULL, ep = 1, approx = FALSE, n = 50, ...){
# 	nmz <- names(x)
# 	
# 	# weights is a variable in x...
# 	if(length(w) == 1){
# 		ii <- NULL
# 		if(is.integer(w)){
# 			ii <- w
# 		}
# 		if(is.character(w)){
# 			ii <- which(names(x)==w)
# 		}
# 		if(!is.null(ii)){
# 			w <- x[,ii]
# 			x <- x[,-ii]
# 		}
# 	}
# 	nd <- ncol(x)
# 	
# 	if(approx){
# 		
# 	cx <- lapply(x, cut, breaks=n,include.lowest=TRUE)
# 	vx <- mapply( function(y,z)  tapply(y,z,mean,nr.rm=TRUE), y = x, z = cx)
# 
# 	M <- matrix(0,nd,nd)
# 	diag(M) <- 1
# 	colnames(M) <- rownames(M) <- nmz
# 	
# 	for(i in 1:(nd-1)){
# 		for(j in (i+1):nd){
# 			z <- subtable(cbind(cx[[i]],cx[[j]]),1:2)
# 			M[i,j] <- M[j, i] <- wdcor(vx[,i][z[,1]],vx[,j][z[,2]],w = z[,3],ep=ep)
# 		}
# 	}
# 	return(M)
# 		
# 	}#approx
# 	
# 	if( any(is.na(x)) ){
# 		x <- na.omit(x)
# 		simpleWarning("NA's found and omitted. Please check.")
# 	}
# 	x <- data.matrix(x)
# 	nd <- ncol(x)
# 	
# 	
# 	if(ncol(x)*nrow(x)^2 > 1e10){
# 	 	simpleWarning("What a big problem. Please use approx = TRUE.")
# 	 	ret <- 42
# 	 	attr(ret,"question") <- "Answer to the Ultimate Question of Life, the Universe, and Everything"
# 	 return(ret)
# 	}
# 	if(is.null(w)){
# 		w <- rep(1,nrow(x))
# 	}
# 	stopifnot(nrow(x) == length(w))
# 	stopifnot(all(w >= 0))
# 		
# 	storage.mode(w) <- "double"
# 	storage.mode(x) <- "double"
# 	
# 	storage.mode(ep) <- "double"
# 	
# 	ret <- .Call("dcorM",x,as.integer(nd),w/sum(w),ep)
# 	M <- matrix(0,nd,nd)
# 	M[lower.tri(M)] <- ret
# 	M <- M + t(M)
#     diag(M) <- 1
# 	colnames(M) <- rownames(M) <- nmz
# 	return(M)
# }


wdcor.data.frame <- function(x, w = NULL, ep = 1, approx = FALSE, n = 50, ...){
  nmz <- names(x)
  # weights is a variable in x...
  if(length(w) == 1){
    ii <- NULL
    if(is.integer(w)){
      ii <- w
    }
    if(is.character(w)){
      ii <- which(names(x)==w)
    }
    if(!is.null(ii)){
      w <- x[,ii]
      x <- x[,-ii]
    }
  }
  nd <- ncol(x)
  
  if(approx){
    cx <- data.table(sapply(x, cut, breaks=n,include.lowest=TRUE, labels = FALSE))
    vx <-
     as.data.frame( mapply(
        function(y,z){ 
          mm <- tapply(y,z,mean,nr.rm=TRUE)
   ret <- rep(0,n)
          ret[as.integer(names(mm))] <- mm
          return(ret)
          },
        y = x, z = cx))
    setkeyv(cx,vars <- names(cx))                 
    if(!is.null(w)){
      cx$w  <- w
      cx <- cx[,list(Freq = sum(w)),by = key(cx)]
    }else{
      cx <- cx[,list(Freq = .N),by = key(cx)]
    }
    M <- diag(1/2,nd)
    dimnames(M) <- list(vars,vars)
    # This to pass the visible binding test... how ugly.
    Freq <- NULL
    M[lower.tri(M)] <- 
      apply(combn(vars,2),2,function(z){
        setkeyv(cx,z)
        cxx <- cx[,list(weight = sum(Freq)),by = key(cx)]
        return(wdcor(
          x = vx[z[1]][cxx[,1,with=FALSE][[1]],],
          y = vx[z[2]][cxx[,2,with=FALSE][[1]],],
          w = cxx[,3,with=FALSE][[1]]))
      })
    M <- M + t(M)
    return(M)
  }#approx
  
  if( any(is.na(x)) ){
    x <- na.omit(x)
    simpleWarning("NA's found and omitted. Please check.")
  }
  x <- data.matrix(x)
  nd <- ncol(x)
  
  
  if(ncol(x)*nrow(x)^2 > 1e10){
    simpleWarning("What a big problem. Please use approx = TRUE.")
    ret <- 42
    attr(ret,"question") <- "Answer to the Ultimate Question of Life, the Universe, and Everything"
    return(ret)
  }
  if(is.null(w)){
    w <- rep(1,nrow(x))
  }
  stopifnot(nrow(x) == length(w))
  stopifnot(all(w >= 0))
  
  storage.mode(w) <- "double"
  storage.mode(x) <- "double"
  
  storage.mode(ep) <- "double"
  
  ret <- .Call("dcorM",x,as.integer(nd),w/sum(w),ep)
  M <- matrix(0,nd,nd)
  M[lower.tri(M)] <- ret
  M <- M + t(M)
  diag(M) <- 1
  colnames(M) <- rownames(M) <- nmz
  return(M)
}