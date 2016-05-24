scoreU.cavgb2 <- function (fac, z, lambda) 
{
	evl <- exp(z %*% lambda)
	ck <- apply(evl, 1, sum) + 1
	pkl <- evl/ck
	pkL <- cbind(pkl, 1 - rowSums(pkl))
	L <- dim(pkL)[2]
	denom <- rowSums(pkL * fac)  # column vector
	num <- fac[, -L]
	midt <- num/as.vector(denom) - 1
	U <- pkl[,-L] * midt
	return(U)
}


scorez.cavgb2 <- function(U,z){
  	n <- dim(U)[1]
  	L1 <- dim(U)[2]  # L-1
  	SC <- matrix(nrow=n)
  	for (l in 1:L1){
   		SCl <-matrix(z*U[,l],nrow= n)  #  U_kl*z_k 
    		colnames(SCl) <- paste(colnames(z),l,sep="")
    		SC <- cbind(SC,SCl)
    	}
  	SC <- SC[,-1]
	
	return(SC)
}

varscore.cavgb2 <- function(SC,w=rep(1,dim(SC)[1])){
    Vsc <- varscore.cgb2(SC,w)
    return(Vsc)
}

desvar.cavgb2 <- function(data=data, SC=SC, ids=NULL, probs=NULL, strata = NULL, 
        variables = NULL, fpc=NULL,
	nest = FALSE, check.strata = !nest, weights=NULL,pps=FALSE,variance=c("HT","YG")) {
	
	desvar.cgb2(data=data, U=SC, ids=ids, probs=probs, strata = strata, variables = variables, 
	fpc=fpc, nest = nest, check.strata = !nest, weights=weights, pps=pps, variance=variance)
}

hess.cavgb2 <- function(U,P,z,w=rep(1, dim(z)[1])){
	dz <- dim(z)	# n x I  	z: matrix of auxiliary variables
	dP <- dim(P)	# n x L  	P: matrix of mixture probabilities
	dU <- dim(U)	# n x (L-1) 	U: matrix of scores, see eq. 21
	if (dP[2] ==2) {
		dU2=1
		dU1=length(U)
	}
	else {
		dU1=dU[1]
		dU2=dU[2]
		}
	L1 <- dP[2]-1
	Lw <- length(w)
	if  ((dU2!=L1)|!((dU1==Lw) | (Lw==1)) | (dz[1]!=dP[1]) | (dz[1]!=dU1) ){
		warning("error in dimensions: no of parameters= ",L1,"; length(w)= ",Lw,"; 
		dim(U)= ",dU[1],",",dU[2],"; dim(P)= ",dP[1],",",dP[2])
		return()
	}
	n  <- dz[1]
	I  <- dz[2]
	L1 <- dU2
	V2 <- matrix(0,nrow=I*L1,ncol=I*L1)
	nn <- expand.grid(gr=colnames(z),par=1:L1)  
	na <- paste(nn$gr,nn$pa,sep="")
	colnames(V2)<- na
	rownames(V2)<- na

	if (L1==1){
		A1 <- -U^2
		A2 <- A1 + U - 2*P[,1]*U
		V2 <- t(z*w)%*%(z*A2)
	}

	else{
	      
		for (i in 1:L1){
		  for (j in 1:L1){
		    a1ij <- -U[,i]*U[,j]
		    a2ij <- -P[,i]*U[,j] - P[,j]*U[,i] + a1ij +(i==j)*U[,i]
		    indi <- (i-1)*I + 1:I
		    indj <- (j-1)*I + 1:I
		    V2[indi,indj] <- t((z*w))%*%(z*a2ij)
		    }
		}
		    
	}
	eigv <- eigen(V2)[[1]]
  	if (max(eigv)>0) {
    		print("Spurious estimates: Fisher information matrix non negative definite.",quote=FALSE)
    		print("Eigenvalues:",quote=FALSE)
    		print(eigv)
  	}
	else{}
	return(V2)
}


vepar.cavgb2 <- function(ml,Vsc, hess)
{
    estimate <- ml[[2]]$par
    V2 <- hess
    V <- solve(V2)
    Vcov <- V %*% Vsc %*% V
    stderr <- sqrt(diag(Vcov))
    Vcor <- diag(1/stderr)%*%Vcov%*%diag(1/stderr)
    names(estimate) <- rownames(Vcov)
    dimnames(Vcor) <- dimnames(Vcov)
    
    return(list(type="parameter",estimate=estimate, stderr=stderr, Vcov=Vcov, Vcor=Vcor))   
}


veind.cavgb2 <- function(group,vepar,shape1, scale, shape2, shape3, pl0, P, decomp="r") {
	L <- length(pl0)
	K <- length(levels(group))
	error <- "FALSE"
	for (k in 1:K){
	  dPk <- length(unique(P[group==levels(group)[k]]))
	  if (dPk>L){
	    error <- "TRUE"
	    warning("the estimated probabilities are not uniquely defined for group ", levels(group)[k])
	  }     
	}
  if (error) return()
	L2 <- L-2
	indic <- list()
	
	for (k in 1:K){
		pk <- as.vector(unique(P[group==levels(group)[k]]))
		esti <- derivind.cgb2(shape1, scale, shape2, shape3, pl0, pk, decomp=decomp) 
		MFDI <- esti[["jacobian"]]
		indi <- c(k+K*(0:L2))
		Vcov.gr <- vepar[["Vcov"]][indi,indi]
		Vcov <- MFDI%*%Vcov.gr%*%t(MFDI)
		stderr <- sqrt(diag(Vcov))
		Vcor <- diag(1/stderr)%*%Vcov%*%diag(1/stderr)
		rownames(Vcor) <- colnames(Vcor) <- rownames(MFDI)
		ngroup <- levels(group)[k]
		indic[[k]] <- list(group=ngroup,estimate=esti[["estimate"]],stderr=stderr,Vcov=Vcov,Vcor=Vcor)
	}
	indic[["type"]] <- "indicator"
    return(indic)
}