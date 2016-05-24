#The function performs biclustering of the data matrix by sparse singular value decomposition with nested stability selection.
#arguments:
#
#	x 				The matrix to be clustered.
#   steps			Number of subsamples used to perform the stability selection
#   pcerv			Per comparsion wise error rate to control the number of falsely selected right singular vector coefficients (columns/samples).
#	pceru			Per comparsion wise error rate to control the number of falsely selected left singular vector coefficients (rows/genes).
#	ss.thr			Range of the cutoff threshold (relative selection frequency) for the stability selection.
#	size			Size of the subsamples used to perform the stability selection.  
#	gamm			Weight parameter for the adaptive LASSO, nonnegative constant (default = 0, LASSO).
#	iter			Maximal number of iterations to fit a single bicluster.
#	nbiclust		Maximal number of biclusters. 
#	merr			Threshold to decide convergence. 
#	cols.nc			Allow for negative correlation of columns (samples) over rows (genes).
#	rows.nc			Allow for negative correlation of rows (genes) over columns (samples).
#	row.overlap		Allow rows to overlap between biclusters. 
#	col.overlap		Allow columns to overlap between biclusters. 
#	row.min			Minimal number of rows.
#	col.min			Minimal number of columns.
#	pointwise		If TRUE performs a fast pointwise stability selection instead of calculating the complete stability path.  
#	start.iter		Number of starting iterations in which the algorithm is not allowed to converge. 
#	savepath		Saves the stability path in order plot the path with the stabpathplot function.
	
#Note that pointwise needs to be TRUE to save the path. For extreme high dimensional data sets (e.g. the lung cancer example) the resulting
#biclust object may exceed the available memory.



s4vd <- function(
		X,
		steps = 100,
		pcerv = 0.1,
		pceru = 0.1,
		ss.thr = c(0.6,0.65),
		size = 0.5,
		gamm = 0,
		iter = 100,
		nbiclust = 10,
		merr = 10^(-3),
		cols.nc=TRUE,
		rows.nc=TRUE,
		row.overlap=TRUE,
		col.overlap=TRUE,
		row.min=1,
		col.min=1,
		pointwise=TRUE,
		start.iter=3,
		savepath=FALSE
){
	MYCALL<-match.call()
	startX <- X
	p.ini <- nrow(X)
	n.ini <- ncol(X)
	rowsin <- rep(TRUE,p.ini)	
	colsin <- rep(TRUE,n.ini)
	stop <- FALSE
	start <- TRUE
	info <- Rows <- Cols <- vc <- uc <- list()
	for(k in 1:nbiclust){
		gc()
		cat("Bicluster",k)
		rows <- rep(FALSE,nrow(startX))
		cols <- rep(FALSE,ncol(startX))
		if(is.null(nrow(X))|is.null(ncol(X))){
			number <- k-1
			stop <- TRUE
			break
		}
		if(nrow(X)==0|ncol(X)==0){
			number <- k-1
			stop <- TRUE
			break
		}
		SVD <- svd(X,nu=1,nv=1)
		v0 <- SVD$v
		u0 <- SVD$u
		d0 <- SVD$d
		vc <- uc <- list()
		if((length(u0)*size)<=2|(length(v0)*size)<=2){
			cat("submatrix to small for resampling","\n")
			number <- k-1
			stop <- TRUE
			break
		}
		if(pointwise){
			for(i in 1:iter){
				if(i > start.iter) start <- FALSE
				uc <- updateu.pw(X,v0,pceru,p.ini,ss.thr,steps,size,gamm,rows.nc,uc$l)
				u1 <- uc[[1]]/sqrt(sum(uc[[1]]^2)) 
				u1[is.na(u1)] <- 0
				vc <- updatev.pw(X,u1,pcerv,n.ini,ss.thr,steps,size,gamm,cols.nc,vc$l)
				v1 <- vc[[1]]/sqrt(sum(vc[[1]]^2)) 
				v1[is.na(v1)] <- 0
				if(uc[[3]] & i > start.iter){
					cat("rows not stable")
					stop <- TRUE
					break}
				if(vc[[3]] & i > start.iter){
					cat("columns not stable")
					stop <- TRUE
					break}
				ud <- sqrt(sum((u0-u1)^2))
				vd <- sqrt(sum((v0-v1)^2))
				#cat("iter: ",i," rows: ",sum(uc$sp>=uc$thr)," cols: ",sum(vc$sp>=uc$thr)
				#		," merr: ",min(c(ud,vd)),"row.thr:",uc[[5]],"col.thr",vc[[5]],"\n")
				cat(".")
				v0 <- v1
				u0 <- u1
				if(min(c(vd,ud)) < merr & i > start.iter)break
			}
		}else{
			for(i in 1:iter){
				if(i > start.iter) start <- FALSE
				uc <- updateu(X,v0,pceru,p.ini,ss.thr,steps,size,gamm,rows.nc,savepath)
				u1 <- uc[[1]]/sqrt(sum(uc[[1]]^2)) 
				u1[is.na(u1)] <- 0
				vc <- updatev(X,u1,pcerv,n.ini,ss.thr,steps,size,gamm,cols.nc,savepath)
				v1 <- vc[[1]]/sqrt(sum(vc[[1]]^2)) 
				v1[is.na(v1)] <- 0
				if(uc[[3]] & i > start.iter){
					cat("rows not stable")
					stop <- TRUE
					break}
				if(vc[[3]] & i > start.iter){
					cat("columns not stable")
					stop <- TRUE
					break}
				ud <- sqrt(sum((u0-u1)^2))
				vd <- sqrt(sum((v0-v1)^2))
				#cat("iter: ",i," rows: ",sum(uc$sp>=uc$thr)," cols: ",sum(vc$sp>=uc$thr)
				#		," merr: ",min(c(ud,vd)),"row.thr:",uc[[5]],"col.thr",vc[[5]],"\n")
				cat(".")
				v0 <- v1
				u0 <- u1
				if(min(c(vd,ud)) < merr & i > start.iter)break
			}
		}
		stableu <- uc$sp >= uc$thr
		stablev <- vc$sp >= vc$thr
		d0 <- as.numeric(t(u0)%*%X%*%v0)
		u0[!stableu] <- 0
		v0[!stablev] <- 0
		rows[rowsin] <- u0!=0
		cols[colsin] <- v0!=0
		Rows[[k]] <- rows
		Cols[[k]] <- cols
		if(stop){
			number <- k-1
			break
		}
		if(i==iter){
			number <- k-1
			stop <- TRUE
			cat("Fail to converge! Increase the number of iterations !","\n")
			gc()
			break
		}
		if(!row.overlap){
			rowsin[rows] <- FALSE
			X <- startX[rowsin,colsin]
			info[[k]] <- list(vc,uc,layer=list(u0,v0,d0))
		} 
		if(!col.overlap){
			colsin[cols] <- FALSE
			X <- startX[rowsin,colsin]
			info[[k]] <- list(vc,uc,layer=list(u0,v0,d0))
		} 
		if(row.overlap&col.overlap){
			temp <- svd(X[rows,cols]) 
			#X <- X - (d0*u0%*%t(v0))
			X[rows,cols] <- X[rows,cols] - (temp$d[1]*temp$u[,1]%*%t(temp$v[,1]))
			info[[k]] <- list(vc,uc,layer=list(u0,v0,d0))
		}
		cat("\n")
	}
	if(!stop) number <- k
	params <- list(steps = steps, pcerv=pcerv, pceru=pceru, iter=iter, ss.thr=ss.thr, size=size, gamm=gamm, row.overlap=row.overlap, col.overlap=col.overlap,
			rows.nc=rows.nc, cols.nc=cols.nc, nbiclust=nbiclust, merr=merr, row.min=row.min, col.min=col.min, pointwise=pointwise, start.iter=start.iter, savepath=savepath, Call=MYCALL)  
	RowxNumber=t(matrix(unlist(Rows),byrow=T,ncol=length(Rows[[1]])))
	NumberxCol=matrix(unlist(Cols),byrow=T,ncol=length(Cols[[1]]))
	if(number)RowxNumber <- matrix(RowxNumber[,1:number],ncol=number)
	if(number)NumberxCol <- matrix(NumberxCol[1:number,],nrow=number)
	Number <- number
	info[[Number+1]] <- params
	return(BiclustResult(params,RowxNumber,NumberxCol,Number,info))
}


#update v
updatev <- function(X,u0,pcer,n.ini,ss.thr,steps,size,gamm,cols.nc=FALSE,savepath=FALSE){
	n.ini  <-	n <- ncol(X)
	err <- pcer*n.ini
	ols <- t(X)%*%u0	
	stop <- FALSE
	lambdas <- sort(c(abs(ols),0),decreasing=TRUE)
	if(savepath) selprobpath <- matrix(nrow=length(ols),ncol=length(lambdas))
	qs <- numeric(length(lambdas))
	thrall <- numeric(length(lambdas))
	ls <- length(lambdas)
	if(cols.nc){
		for(l in 1:length(lambdas)){
			temp <- adaLasso.nc(t(X),u0,lambdas[l],steps,size,gamm)
			t <- temp!=0
			qs[l] <- mean(colSums(t))
			sp <- rowMeans(t)
			if(savepath)selprobpath[,l] <- sp
			thrall[l] <- ((qs[l]^2/(err*n.ini))+1)/2
			if(thrall[l]>=ss.thr[1]){
				ls <- l
				break
			}
		}
	}else{
		for(l in 1:length(lambdas)){
			temp <- adaLasso(t(X),u0,lambdas[l],steps,size,gamm)
			t <- temp!=0
			qs[l] <- mean(colSums(t))
			sp <- rowMeans(t)
			if(savepath)selprobpath[,l] <- sp
			thrall[l] <- ((qs[l]^2/(err*n.ini))+1)/2
			if(thrall[l]>=ss.thr[1]){
				ls <- l
				break
			}
		}
	}
	thr <- thrall[ls]
	if(thr > ss.thr[2]){
		while(pcer <= 0.5){
			pcer <- pcer + 0.01	
			thrall <- ((qs^2/((pcer*n.ini)*n.ini))+1)/2
			thr <- thrall[ls]
			if(thr < ss.thr[2])break
		}
	} 
	stable <- which(sp>=thr)
	if(length(stable)==0)stop <- TRUE
	vc <- rep(0,n)
	delta <- lambdas[ls]/(abs(ols)^gamm)  
	vc <- (sign(ols)*(abs(ols)>=delta)*(abs(ols)-delta))
	if(savepath){
		return(list(vc=vc,sp=sp,stop=stop,qs=qs,thr=thr,l=ls,delta=delta,selprobpath=selprobpath))
	}else{
		return(list(vc=vc,sp=sp,stop=stop,qs=qs,thr=thr,l=ls,delta=delta))
	}
}

#update u
updateu <- function(X,v0,pcer,p.ini,ss.thr,steps,size,gamm,rows.nc=FALSE,savepath=FALSE,start=FALSE){
	p.ini <- p <- nrow(X)
	err <- pcer*p.ini
	ols <- X%*%v0
	stop <- FALSE
	lambdas <- sort(c(abs(ols),0),decreasing=TRUE)
	if(savepath) selprobpath <- matrix(nrow=length(ols),ncol=length(lambdas))
	qs <- numeric(length(lambdas))
	thrall <- numeric(length(lambdas))
	ls <- length(lambdas)
	if(rows.nc){
		for(l in 1:length(lambdas)){
			temp <- adaLasso.nc(X,v0,lambdas[l],steps,size,gamm)
			t <- temp!=0
			qs[l] <- mean(colSums(t))
			sp <- rowMeans(t)
			if(savepath)selprobpath[,l] <- sp
			thrall[l] <- ((qs[l]^2/(err*p.ini))+1)/2
			if(thrall[l]>=ss.thr[1]){
				ls <- l
				break
			}
		}
	}else{
		for(l in 1:length(lambdas)){
			temp <- adaLasso(X,v0,lambdas[l],steps,size,gamm)
			t <- temp!=0
			qs[l] <- mean(colSums(t))
			sp <- rowMeans(t)
			if(savepath)selprobpath[,l] <- sp
			thrall[l] <- ((qs[l]^2/(err*p.ini))+1)/2
			if(thrall[l]>=ss.thr[1]){
				ls <- l
				break
			}
		}
	}
	thr <- thrall[ls]
	if(thr > ss.thr[2]){
		while(pcer <= 0.5){
			pcer <- pcer + 0.01	
			thrall <- ((qs^2/((pcer*p.ini)*p.ini))+1)/2
			thr <- thrall[ls]
			if(thr < ss.thr[2])break
		}
	} 
	stable <- which(sp>=thr)
	if(length(stable)==0)stop <- TRUE
	uc <- rep(0,p)
	delta <- lambdas[ls]/(abs(ols)^gamm)  
	uc <- (sign(ols)*(abs(ols)>=delta)*(abs(ols)-delta))
	if(savepath){
		return(list(uc=uc,sp=sp,stop=stop,qs=qs,thr=thr,l=ls,delta=delta,selprobpath=selprobpath))
	}
	else{
		return(list(uc=uc,sp=sp,stop=stop,qs=qs,thr=thr,l=ls,delta=delta))
	}	
}

#update v pointwise

updatev.pw <- function(X,u0,pcer,n.ini,ss.thr,steps,size,gamm,cols.nc=FALSE,l=NULL,start=FALSE){
	n.ini <- n <- ncol(X)
	err <- pcer*n.ini
	ols <- t(X)%*%u0
	stop <- FALSE
	lambdas <- sort(c(abs(ols),0),decreasing=TRUE)
	qs <- numeric(length(lambdas))
	thrall <- numeric(length(lambdas))
	ls <- length(lambdas)
	if(is.null(l)) l <- which(lambdas==quantile(lambdas,0.5,type=1))[1]
	#search for a lambda
	l.min <- 1
	l.max <- length(lambdas)
	if(cols.nc){
		for(g in 1:(length(lambdas))){
			temp <- adaLasso.nc(t(X),u0,lambdas[l],steps,size,gamm)
			t <- temp!=0
			qs[l] <- mean(colSums(t))
			sp <- rowMeans(t)
			thrall[l] <- ((qs[l]^2/(err*n.ini))+1)/2
			if(thrall[l] >= ss.thr[1] & thrall[l] <= ss.thr[2]){
				ls <- l
				break
			} 
			if(thrall[l] < ss.thr[1]){          
				l.min <- l
				if(l == length(lambdas))break   
				if(thrall[l+1]> ss.thr[2]){
					ls <- l+1
					if(thrall[l+1]>1) ls <-l
					temp <- adaLasso.nc(t(X),u0,lambdas[ls],steps,size,gamm)
					t <- temp!=0
					qs[ls] <- mean(colSums(t))
					sp <- rowMeans(t)
					thrall[ls] <- ((qs[ls]^2/(err*n.ini))+1)/2
					break
				}
				l <- min(length(lambdas),l.max,l + ceiling(length(lambdas)/(g+1))) 
				while(thrall[l]!=0){  
					l <- l-1
					if(l == 0)break
				} 
			}
			if(thrall[l] > ss.thr[2]){ 
				l.max <- l
				if(l == 0)break
				if(thrall[l-1]<ss.thr[1]&thrall[l-1]!=0){
					ls <- l
					break
				}
				l <- max(1,l.min,l - ceiling(length(lambdas)/(g+1)))
				while(thrall[l]!=0 ){ 
					l <- l+1
					if(l == length(l))break
				} 
			}
		}
	}else{
		for(g in 1:(length(lambdas))){
			temp <- adaLasso(t(X),u0,lambdas[l],steps,size,gamm)
			t <- temp!=0
			qs[l] <- mean(colSums(t))
			sp <- rowMeans(t)
			thrall[l] <- ((qs[l]^2/(err*n.ini))+1)/2
			if(thrall[l] >= ss.thr[1] & thrall[l] <= ss.thr[2]){
				ls <- l
				break
			} 
			if(thrall[l] < ss.thr[1]){          
				l.min <- l
				if(l == length(lambdas))break   
				if(thrall[l+1]> ss.thr[2]){
					ls <- l +1 
					if(thrall[l+1]>1) ls <-l
					temp <- adaLasso(t(X),u0,lambdas[ls],steps,size,gamm)
					t <- temp!=0
					qs[ls] <- mean(colSums(t))
					sp <- rowMeans(t)
					thrall[ls] <- ((qs[ls]^2/(err*n.ini))+1)/2
					break
				}
				l <- min(length(lambdas),l.max,l + ceiling(length(lambdas)/(g+1))) 
				while(thrall[l]!=0 ){  
					l <- l-1
					if(l == 0)break
				} 
			}
			if(thrall[l] > ss.thr[2]){ 
				l.max <- l
				if(l == 0)break
				if(thrall[l-1]<ss.thr[1]&thrall[l-1]!=0){
					ls <- l
					break
				}
				l <- max(1,l.min,l - ceiling(length(lambdas)/(g+1)))
				while(thrall[l]!=0 ){ 
					l <- l+1
					if(l == length(lambdas))break
				} 
			}
		}
	}	
	#thr <- thrall[ls]
	#if(thr > ss.thr[2]){
	#	while(pcer <= 0.5){
	#		pcer <- pcer + 0.01	
	#		thrall <- ((qs^2/((pcer*n.ini)*n.ini))+1)/2
	#		thr <- thrall[ls]
	#		if(thr < ss.thr[2])	break
	#	}
	#}
	thr <- ((qs[ls]^2/((pcer*n.ini)*n.ini))+1)/2
	stable <- which(sp>=thr)
	if(length(stable)==0)stop <- TRUE
	vc <- rep(0,n)
	delta <- lambdas[l]/(abs(ols)^gamm)  
	vc <- (sign(ols)*(abs(ols)>=delta)*(abs(ols)-delta))
	return(list(vc=vc,sp=sp,stop=stop,qs=qs,thr=thr,l=ls,delta=delta))
}


#update u pointwise
updateu.pw <- function(X,v0,pcer,p.ini,ss.thr,steps,size,gamm,rows.nc=FALSE,l=NULL,start=FALSE){
	p.ini <- p <- nrow(X)
	err <- pcer*p.ini
	ols <- X%*%v0
	stop <- FALSE
	lambdas <- sort(c(abs(ols),0),decreasing=TRUE)
	qs <- numeric(length(lambdas)) 
	thrall <- numeric(length(lambdas))
	ls <- length(lambdas)
	if(is.null(l)) l <- which(lambdas==quantile(lambdas,0.5,type=1))[1]
	#search for a lambda
	l.min <- 1
	l.max <- length(lambdas)
	if(rows.nc){
		for(g in 1:(length(lambdas))){
			temp <- adaLasso.nc(X,v0,lambdas[l],steps,size,gamm)
			t <- temp!=0
			qs[l] <- mean(colSums(t))
			sp <- rowMeans(t)
			thrall[l] <- ((qs[l]^2/(err*p.ini))+1)/2
			if(thrall[l] >= ss.thr[1] & thrall[l] <= ss.thr[2]){
				ls <- l
				break
			} 
			if(thrall[l] < ss.thr[1]){
				l.min <- l
				if(l == length(lambdas))break
				if(thrall[l+1]> ss.thr[2]){
					ls <- l+1
					if(thrall[l+1]>1) ls <-l
					temp <- adaLasso.nc(X,v0,lambdas[ls],steps,size,gamm)
					t <- temp!=0
					qs[ls] <- mean(colSums(t))
					sp <- rowMeans(t)
					thrall[ls] <- ((qs[ls]^2/(err*p.ini))+1)/2
					break
				}
				l <- min(length(lambdas),l.max,l + ceiling(length(lambdas)/(g+1)))
				while(thrall[l]!=0 ){  # if thr for current lambda available decrease lambda 
					l <- l-1
					if(l == 0)break
				} 
			}
			if(thrall[l] > ss.thr[2]){ 
				l.max <- l
				if(l == 0)break
				if(thrall[l-1]<ss.thr[1]&thrall[l-1]!=0){
					ls <- l
					break
				}
				l <- max(1,l.min,l - ceiling(length(lambdas)/(g+1)))
				while(thrall[l]!=0){ 
					l <- l+1
					if(l == length(l))break
				} 
			}
		}
	}else{
		for(g in 1:(length(lambdas))){
			temp <- adaLasso(X,v0,lambdas[l],steps,size,gamm)
			t <- temp!=0
			qs[l] <- mean(colSums(t))
			sp <- rowMeans(t)
			thrall[l] <- ((qs[l]^2/(err*p.ini))+1)/2
			if(thrall[l] >= ss.thr[1] & thrall[l] <= ss.thr[2]){
				ls <- l
				break
			} 
			if(thrall[l] < ss.thr[1]){
				l.min <- l
				if(l == length(lambdas))break
				if(thrall[l+1]> ss.thr[2]){
					ls <- l+1
					if(thrall[l+1]>1) ls <-l
					temp <- adaLasso(X,v0,lambdas[ls],steps,size,gamm)
					t <- temp!=0
					qs[ls] <- mean(colSums(t))
					sp <- rowMeans(t)
					thrall[ls] <- ((qs[ls]^2/(err*p.ini))+1)/2
					break
				}
				l <- min(length(lambdas),l.max,l + ceiling(length(lambdas)/(g+1)))
				while(thrall[l]!=0 ){  
					l <- l-1
					if(l == 0)break
				} 
			}
			if(thrall[l] > ss.thr[2]){ 
				l.max <- l
				if(l == 0)break
				if(thrall[l-1]<ss.thr[1]&thrall[l-1]!=0){
					ls <- l
					break
				}
				l <- max(1,l.min,l - ceiling(length(lambdas)/(g+1)))
				while(thrall[l]!=0){ 
					l <- l+1
					if(l == length(lambdas))break
				} 
			}
		}
	}
	#thr <- thrall[l]
	#if(thr > ss.thr[2]){
	#	while(pcer <= 0.5){
	#		pcer <- pcer + 0.01	
	#		thrall <- ((qs^2/((pcer*p.ini)*p.ini))+1)/2
	#		thr <- thrall[ls]
	#		if(thr < ss.thr[2])break
	#	}
	#}
	thr <- ((qs[ls]^2/((pcer*p.ini)*p.ini))+1)/2
	stable <- which(sp>=thr)
	if(length(stable)==0)stop <- TRUE
	uc <- rep(0,p)
	delta <- lambdas[l]/(abs(ols)^gamm)  
	uc <- (sign(ols)*(abs(ols)>=delta)*(abs(ols)-delta))
	return(list(uc=uc,sp=sp,stop=stop,qs=qs,thr=thr,l=ls,delta=delta))
}

#adaptive Lasso 
adaLasso.nc <- function(X,b,lambda,steps,size,gamm=0){
	subsets <- sapply(1:steps,function(x){sample(1:length(b),length(b)*size)})
	res <- sapply(1:steps,adaLassoSteps.nc,subsets,X,b,lambda,gamm)
	return(res)
}
adaLassoSteps.nc <- function(index,subsets,X,b,lambda,gamm){
	ols <- X[,subsets[,index]]%*%b[subsets[,index]]
	delta <- lambda/(abs(ols)^gamm)                        
	ols <- sign(ols)*(abs(ols)>=delta)*(abs(ols)-delta)
	ols[is.na(ols)] <- 0
	return(ols)
}
adaLasso <- function(X,b,lambda,steps,size,gamm=0){
	subsets <- sapply(1:steps,function(x){sample(1:length(b),length(b)*size)})
	res <- sapply(1:steps,adaLassoSteps,subsets,X,b,lambda,gamm)
	return(res)
}
adaLassoSteps <- function(index,subsets,X,b,lambda,gamm){
	ols <- X[,subsets[,index]]%*%b[subsets[,index]]
	delta <- lambda/(abs(ols)^gamm)                        
	ols <- sign(ols)*(abs(ols)>=delta)*(abs(ols)-delta)
	ols[is.na(ols)] <- 0
	mostof <- sign(sum(sign(ols)))
	if(mostof==0) mostof <- 1
	ols[which(sign(ols) != mostof) ] <- 0
	ols[is.na(ols)] <- 0
	return(ols)
}
