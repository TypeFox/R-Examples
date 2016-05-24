A.mat <- function(X,min.MAF=NULL,max.missing=NULL,impute.method="mean",tol=0.02,n.core=1,shrink=FALSE,return.imputed=FALSE){

impute.EM <- function(W, cov.mat, mean.vec) {
	n <- nrow(W)
	m <- ncol(W)
	S <- matrix(0,n,n)
	for (i in 1:m) {
		Wi <- matrix(W[,i],n,1)
		missing <- which(is.na(Wi))
		if (length(missing) > 0) {
			not.NA <- setdiff(1:n,missing)
			Bt <- solve(cov.mat[not.NA,not.NA],cov.mat[not.NA,missing])
			Wi[missing] <- mean.vec[missing] + crossprod(Bt,Wi[not.NA]-mean.vec[not.NA])
			C <- cov.mat[missing,missing] - crossprod(cov.mat[not.NA,missing],Bt)
			D <- tcrossprod(Wi)
			D[missing,missing] <- D[missing,missing] + C
			W[,i] <- Wi
		} else {D <- tcrossprod(Wi)}
		S <- S + D
	}	
	return(list(S=S,W.imp=W))
}

cov.W.shrink <- function(W) {
	m <- ncol(W)
	n <- nrow(W)
	Z <- t(scale(t(W),scale=FALSE))
	Z2 <- Z^2
	S <- tcrossprod(Z)/m
	target <- mean(diag(S))*diag(n)
	var.S <- tcrossprod(Z2)/m^2-S^2/m
	b2 <- sum(var.S)
	d2 <- sum((S-target)^2)
	delta <- max(0,min(1,b2/d2))
	print(paste("Shrinkage intensity:",round(delta,2)))
	return(target*delta + (1-delta)*S)
}

X <- as.matrix(X)
n <- nrow(X)
frac.missing <- apply(X,2,function(x){length(which(is.na(x)))/n})
missing <- max(frac.missing) > 0
freq <- apply(X + 1, 2, function(x) {mean(x, na.rm = missing)})/2
MAF <- apply(rbind(freq,1-freq),2,min)
if (is.null(min.MAF)) {min.MAF <- 1/(2*n)}
if (is.null(max.missing)) {max.missing <- 1 - 1/(2*n)}
markers <- which((MAF >= min.MAF)&(frac.missing <= max.missing)) 
m <- length(markers)
var.A <- 2 * mean(freq[markers] * (1 - freq[markers]))
one <- matrix(1, n, 1)

mono <- which(freq*(1-freq)==0)
X[,mono] <- 2*tcrossprod(one,matrix(freq[mono],length(mono),1))-1

freq.mat <- tcrossprod(one, matrix(freq[markers], m, 1))
W <- X[, markers] + 1 - 2 *freq.mat 

if (!missing) {
    if (shrink) {
		W.mean <- rowMeans(W)
		cov.W <- cov.W.shrink(W)
		A <- (cov.W+tcrossprod(W.mean))/var.A	
	} else {
		A <- tcrossprod(W)/var.A/m	
	}
	rownames(A) <- rownames(X)
	colnames(A) <- rownames(A)
	if (return.imputed) {
		return(list(A=A,imputed=X))		
	} else {
		return(A)
	}
} else {
    #impute
    isna <- which(is.na(W)) 
	W[isna] <- 0
	
	if (toupper(impute.method)=="EM") {
		if (m < n) {
			print("Linear dependency among the lines: imputing with mean instead of EM algorithm.")
		} else {
			mean.vec.new <- matrix(rowMeans(W),n,1)
			cov.mat.new <- cov(t(W))
			if (qr(cov.mat.new)$rank < nrow(cov.mat.new)-1) {
				print("Linear dependency among the lines: imputing with mean instead of EM algorithm.")
			} else {

			#do EM algorithm
		    W[isna] <- NA
			A.new <- (cov.mat.new + tcrossprod(mean.vec.new))/var.A
			err <- tol+1
			print("A.mat converging:")
			while (err >= tol) {
				A.old <- A.new
				cov.mat.old <- cov.mat.new
				mean.vec.old <- mean.vec.new
				if ((n.core > 1) & requireNamespace("parallel",quietly=TRUE)) {
 					it <- split(1:m,factor(cut(1:m,n.core,labels=FALSE)))
					pieces <- parallel::mclapply(it,function(mark2){impute.EM(W[,mark2],cov.mat.old,mean.vec.old)},mc.cores=n.core)
				} else {
					pieces <- list()
					pieces[[1]] <- impute.EM(W,cov.mat.old,mean.vec.old)
				}
				n.pieces <- length(pieces)
				S <- matrix(0,n,n)
				W.imp <- numeric(0)
				for (i in 1:n.pieces) {
					S <- S + pieces[[i]]$S
					W.imp <- cbind(W.imp,pieces[[i]]$W.imp)
				}
				mean.vec.new <- matrix(rowMeans(W.imp),n,1)
				cov.mat.new <- (S-tcrossprod(mean.vec.new)*m)/(m-1)
				A.new <- (cov.mat.new + tcrossprod(mean.vec.new))/var.A
				err <- norm(A.old-A.new,type="F")/n
				print(err,digits=3)
			}
			rownames(A.new) <- rownames(X)
			colnames(A.new) <- rownames(A.new)

			if (return.imputed) {
				Ximp <- W.imp - 1 + 2*freq.mat
				colnames(Ximp) <- colnames(X)[markers]
				rownames(Ximp) <- rownames(X)
				return(list(A=A.new,imputed=Ximp))
			} else {
				return(A.new)
			}
		  	} #else EM 
  		} #else EM
  	} #else EM
  		
  	#imputing with mean
	if (shrink) {
		W.mean <- rowMeans(W)
		cov.W <- cov.W.shrink(W)
		A <- (cov.W+tcrossprod(W.mean))/var.A	
	} else {
		A <- tcrossprod(W)/var.A/m	
	}
	rownames(A) <- rownames(X)
	colnames(A) <- rownames(A)

	if (return.imputed) {
		Ximp <- W - 1 + 2*freq.mat
		colnames(Ximp) <- colnames(X)[markers]
		rownames(Ximp) <- rownames(X)
		return(list(A=A,imputed=Ximp))		
	} else {
		return(A)
	}
} #else missing

} #A.mat

