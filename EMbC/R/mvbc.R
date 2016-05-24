# The EMbC Package for R
#
# Copyright 2013, 2014, 2015 Joan Garriga <jgarriga@ceab.csic.es>, Aitana Oltra <aoltra@ceab.csic.es>, John R.B. Palmer <johnrbpalmer@gmail.com>, Frederic Bartumeus <fbartu@ceab.csic.es>
#
#   EMbC is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or  (at your option) any later version.
#
# EMbC is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with this program.  If not, see http://www.gnu.org/licenses.



# bivariate EM behavioral Clustering Algorithm
# --------------------------------------------

setGeneric("clst",function(bC,maxItr=200,info=0,check=0,...){standardGeneric("clst")})

setMethod("clst",signature("binClst"),function(bC,maxItr,info,...){

	bC@m <- ncol(bC@X)
	bC@k <- 2**bC@m
	bC@n <- nrow(bC@X)
	bC@C <- getColors(bC@k)

	minVar <- bC@stdv**2
	lblDst <- numeric()

	# delimiters
	Rdel <- lapply(1:bC@k,function(k){
			lapply(1:bC@m,function(m){
				bk <- as.integer(rev(intToBits(k-1)[1:bC@m]))
				if (bk[m]==0) {
					kL <- strtoi(paste(bk,collapse=''),base=2)+1
					bk[m] <- 1
					kH <- strtoi(paste(bk,collapse=''),base=2)+1
					list(m,kL,kH)}
				})
			})
	Rdel <- t(matrix(unlist(Rdel),c(3,bC@k)))

	bC@R <- getRmatrix(bC,info)
	rownames(bC@R) <- getkLbls(bC,kNmbrs=TRUE)
	colnames(bC@R) <- c(sapply(seq(bC@m),function(m) paste('X',m,'.min',sep='')),sapply(seq(bC@m),function(m) paste('X',m,'.max',sep='')))

	# initial assignements
	bC@W <- matrix(rep(1/bC@k,bC@n*bC@k),c(bC@n,bC@k))
	if (maxItr==0) bC@A <- getClusters(bC)
	else bC@A <- floor(runif(bC@n,min=1,max=bC@k+1))

	# parameters
	bC@P <- lapply(seq(bC@k),function(k) list(M=rep(0,bC@m),S=matrix(rep(0,bC@m+bC@m**2),c(bC@m,bC@m+1))))
	bC@P <- getTheta(bC,minVar)

	# show initial status
	stepInfo(bC,-0,length(bC@A),0,info)

	# optimization loop
	A_ <- bC@A
	R_ <- bC@R
	stableItr <- 0

	if (maxItr==0) return(bC)
	for (itrNmb in 1:maxItr){

		# likelihood densities
		L <- array(-Inf,c(bC@n,bC@k))
		for (k in seq(bC@k)){
			if (mean(bC@W[,k])>0){
				Xb <- bound2R(bC,k)
				if (length(which(Xb==TRUE))>1){
					UW  <- bC@U*matrix(rep(bC@W[,k],bC@m),c(bC@n,bC@m))
					bC@P[[k]] <- gssTheta(bC,Xb,UW,k,minVar)
					if (!is.null(bC@P[[k]]$S))
						L[,k]  <- log(mean(bC@W[,k])) + dmnorm(bC@X,bC@P[[k]]$M,bC@P[[k]]$S,log=TRUE)
					}
				}
			}

		# densities2weights
		bC@W <- dens2wght(L)

		# new clustering assignments
		bC@A <- getClusters(bC)

		# likelihood
    	bC@L[itrNmb] <- getLkh(bC,L)

		# show current status
		if (info>=0) stepInfo(bC,bC@L[length(bC@L)],length(which(bC@A!=A_)),itrNmb,info)

		# check Cycles and stop criterium
		if (itrNmb>1){
			lblDst[itrNmb] <- sqrt(sum((bC@A-A_)^2))
			if (lblDst[itrNmb]==0 && all(bC@R==R_)){
				stableItr <- stableItr +1
				if (stableItr==2){
					print('... Stable clustering',quote=FALSE)
					break
					}
				}
			else stableItr <- 0
			}
		if (itrNmb>100){
			if (chckCycl(lblDst)){
				print('... break: optimization cycle',quote=FALSE)
				break
				}
			}

		A_ <- bC@A
		R_ <- bC@R

		# M-step: New Reference Values
		for (i in 1:dim(Rdel)[1]){
			m  <- Rdel[i,1]
			kL <- Rdel[i,2]
			kH <- Rdel[i,3]
			r <- newDel(bC,m,kL,kH)
			if (!is.null(r)){
				bC@R[kL,bC@m+m] <- r
				bC@R[kH,m] <- r}
			}

		}
	return(bC)
	})

# EMbC internal functions
# -----------------------

# Initial Reference Values
getRmatrix <- function(bC,info){

  if (info>=0) cat('... computing starting delimiters:\n')

  R <- matrix(rep(0,bC@k*bC@m*2),c(bC@k,bC@m*2))
  remainingSplits <- list(list(vars=1:bC@m,data=bC@X,clst=rep(-1,bC@m)))

  while (length(remainingSplits)>0){

    remainingVars <- remainingSplits[[1]]$vars
    splittingData <- remainingSplits[[1]]$data

	hlf <- ceiling(dim(splittingData)[1]/2)
	dLst <- lapply(remainingVars,function(m){
		xLst <- as.numeric(names(table(splittingData[,m])))
		nLst <- sapply(xLst, function(s) length(which(splittingData[,m]<=s)))
		s <- xLst[which.min(abs(nLst-hlf))]
		return(list(m=m,splitVal=s,kPrc=sum((splittingData[,m]-s)**2)/sum(abs(splittingData[,m]-s))**2))})

	splitIdx <-  which.max(lapply(dLst,function(d) d$kPrc))
    splitVar <- dLst[[splitIdx]]$m
    splitVal <- dLst[[splitIdx]]$splitVal

    clusterKey <- remainingSplits[[1]]$clst
    patternKey <- which(clusterKey!=-1)
    for (k in 1:bC@k){
      bk <- as.integer(rev(intToBits(k-1)[1:bC@m]))
      if (all(bk[patternKey]==clusterKey[patternKey])){
        if (bk[splitVar]==0){
          R[k,splitVar] <- min(bC@X[,splitVar])
          R[k,(splitVar+bC@m)] <- splitVal
		  if (info>=0){
			  del <- getkLbls(bC)[[k]]
			  substr(del,splitVar,splitVar) <- 'l'
			  cat('... delimiter ',del,': ',splitVal,'\n')
			  }
		  }
        else{
          R[k,splitVar] <- splitVal
          R[k,(splitVar+bC@m)] <- max(bC@X[,splitVar])
		  if (info>=0){
			  del <- getkLbls(bC)[[k]]
			  substr(del,splitVar,splitVar) <- 'h'
			  cat('... delimiter ',del,': ',splitVal,'\n')
			  }
		  }
        }
      }

    remainingVars <- remainingVars[remainingVars!=splitVar]
    if (length(remainingVars)>0){
	  lowData <- which(splittingData[,splitVar]<=splitVal)
      if (length(lowData)>1){
		lowClst <- clusterKey
		lowClst[splitVar] <- 0
		remainingSplits[[length(remainingSplits)+1]] <- list(vars=remainingVars,data=splittingData[lowData,],clst=lowClst)
		}
      hghData <- which(splittingData[,splitVar]>=splitVal)
      if (length(hghData)>1){
        hghClst <- clusterKey
        hghClst[splitVar] <- 1
		remainingSplits[[length(remainingSplits)+1]] <- list(vars=remainingVars,data=splittingData[hghData,],clst=hghClst)
		}
      }

    remainingSplits[[1]] <- NULL
    }
  return(R)
  }

# cluster assignments
getClusters <- function(bC) {
	A <- lapply(1:bC@n,function(i){
		if (any(bC@U[i,]==0)) c <- (bC@k+1)
		else{
			c <- which(bC@W[i,]==max(bC@W[i,]))
			if (length(c)==1) c
			else {
				K <- lapply(c,function(k) {all(bC@X[i,]>=bC@R[k,1:bC@m])&all(bC@X[i,]<=bC@R[k,(bC@m+1):(2*bC@m)])})
				if (all(K==FALSE)||length(which(K==TRUE))>1) bC@k+1
				else c[which(K==TRUE)]
				}
			}
		})
	return(as.numeric(A))
	}

# Rbounded sets
bound2R <- function(bC,k) {
	return(apply(bC@X,1,function(x){all(x>=bC@R[k,1:bC@m])&all(x<=bC@R[k,(bC@m+1):(2*bC@m)])}))
	}

# (Rbounded) gaussian parameters
gssTheta <- function(bC,Xb,UW,k,minVar,lower=-Inf,upper=Inf) {
	Xk  <- bC@X[which(Xb==TRUE),]
	UWk <- UW[which(Xb==TRUE),]
	Mk <- as.numeric(lapply(1:bC@m,function(m){sum(UWk[,m]*Xk[,m])/sum(UWk[,m])}))
	Sk <- matrix(rep(0,bC@m**2),c(bC@m,bC@m))
	for (mi in 1:bC@m){
		if (!is.na(Mk[mi])){
			Sk[mi,mi:bC@m] <- sapply(mi:bC@m,function(mj){
#~ 				UWij <- sqrt(UWk[,mi]**2+UWk[,mj]**2)/sqrt(2)
#~ 				ifelse ((sum(UWij)<=0),-Inf,sum(UWij*(Xk[,mi]-Mk[mi])*(Xk[,mj]-Mk[mj]))/sum(UWij))
				UWij <- sqrt(UW[,mi]**2+UW[,mj]**2)/sqrt(2)
				ifelse ((sum(UWij)<=0),-Inf,sum(UWij*(bC@X[,mi]-Mk[mi])*(bC@X[,mj]-Mk[mj]))/sum(UWij))
				})
			}
		}
	for (mi in 1:bC@m){
		if (Sk[mi,mi]>0 && Sk[mi,mi]<minVar[mi]){
			Sk[mi,] <- rep(0,bC@m)
			Sk[mi,mi] <- minVar[mi]
			}
		Sk[,mi] <- Sk[mi,]
		}
	if (is.null(pd.solve(Sk,silent=TRUE))){
		cat('... singular covariance matrix for',sprintf("cluster=%d",k),'\n')
		Sk <- NULL
		}
	return(list(M=Mk,S=Sk))
	}

# densities to weights conversion
dens2wght <- function(L){
	return(t(apply(L,1,function(x) exp(x-max(x))/sum(exp(x-max(x))))))
	}

# likelihood computation
getLkh  <- function(bC,L){
	Lmax <- apply(L,1,function(l) max(l))
	L <- t(apply(cbind(L,Lmax),1,function(l) exp(l-l[length(l)])))[,1:ncol(L)]
	return(sum(Lmax+log(apply(L,1,sum)))/bC@n)
	}

# delimiters
newDel <- function(bC,m,kL,kH){
	subSet <- which(bC@X[,m]>=bC@P[[kL]]$M[m] & bC@X[,m]<=bC@P[[kH]]$M[m])
	if (mean(bC@W[,kL])==0 || is.null(bC@P[[kL]]$S) || length(subSet)<2)
		if (length(which(bC@A==kL))==0)
			r <- min(bC@X[,m])
		else
			r <- bC@P[[kL]]$M[m]
	else if (mean(bC@W[,kH])==0 || is.null(bC@P[[kH]]$S) || length(subSet)<2)
		if (length(which(bC@A==kH))==0)
			r <- max(bC@X[,m])
		else
			r <- bC@P[[kH]]$M[m]
	else {
		X <- bC@X[subSet,]
		for (j in seq(bC@m))
			if (j!=m) X[,j] <- bC@P[[kL]]$M[j] + (X[,m]-bC@P[[kL]]$M[m])*(bC@P[[kH]]$M[j]-bC@P[[kL]]$M[j])/(bC@P[[kH]]$M[m]-bC@P[[kL]]$M[m])
		# cluster lognormal densities
		dL <- log(mean(bC@W[,kL])) + dmnorm(X,bC@P[[kL]]$M,bC@P[[kL]]$S,log=TRUE)
		dH <- log(mean(bC@W[,kH])) + dmnorm(X,bC@P[[kH]]$M,bC@P[[kH]]$S,log=TRUE)
		# densities to weights
		W <- dens2wght(cbind(dL,dH))
		# find equiprobability point
		if (all(W[,1]==W[,2]))
			return(NULL)
		else if (all(W[,1]<W[,2]))
			if (length(which(bC@A==kL))==0)
				r <- min(bC@X[,m])
			else
				r <- bC@P[[kL]]$M[m]
		else if (all(W[,1]>W[,2]))
			if (length(which(bC@A==kH))==0)
				r <- max(bC@X[,m])
			else
				r <- bC@P[[kH]]$M[m]
		else{
			r <- X[which.min(abs(W[,1]-W[,2])),m]
			r <- max(r,bC@P[[kL]]$M[m])
			r <- min(r,bC@P[[kH]]$M[m])
			}
		return(r)
		}
	}

# show current status
stepInfo <- function(bC,lkh,lblDff,itrNmb,info) {
	it <- formatC(itrNmb,width=3)
	lk <- formatC(lkh,width=10,format='e',flag="-")
	df <- formatC(lblDff,width=8)
	nK <- formatC(sum(as.numeric(lapply(1:bC@k,function(k) as.integer(any(bC@A==k))))),width=6)
	A <- as.numeric(lapply(1:bC@k,function(k){formatC(round(length(bC@A[bC@A==k])/bC@n,4),width=6,flag="-")}))
	print(paste(it,lk,nK,df,sep="  "),quote=FALSE)
   	if (info>0) stts(bC)
   	if (info>1)	print(bC@R)
	}

# check cycles
chckCycl <- function(lblDst){
	for (i in 1:10){
		cue <- seq((length(lblDst)-i),length(lblDst))
		cueMean <- mean(lblDst[cue])
		cueSmDv <- sum(lblDst[cue]-cueMean)
		if (all(unlist(lapply(1:5,function(j){(sum(lblDst[(cue-((j*i)+1))]-cueMean)==cueSmDv)}))))
			return(TRUE)
		}
	return(FALSE)
	}

# compute parameters
# (also intended for external methods, i.e. "rlbl", allthough not used at the moment)
getTheta <- function(bC,minVar){
	for (k in seq(bC@k)){
		if (mean(bC@W[,k])>0){
			Xb <- bound2R(bC,k)
			if (length(which(Xb==TRUE))>1){
				UW  <- bC@U*matrix(rep(bC@W[,k],bC@m),c(bC@n,bC@m))
				bC@P[[k]] <- gssTheta(bC,Xb,UW,k,minVar)
				}
			}
		}
	return(bC@P)
	}
