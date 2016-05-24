
###extract the first fw PCs such that they explain at least cutoff*100%
### of total variation from X;
### Input: X: nObs by nSNPs;
###        cutoff: threshold to determine how many first few PCs to use;
### Output: a matrix consisting of the first few PCs;

extractPCs<-function(X, cutoff=0.95){

    if (is.null(ncol(X)) || ncol(X)<2) Xpcs<-X
    else{
        X.pca<-prcomp(X, center=FALSE)
        ##find the min k s.t. the first k PCs explain >= cutoff variations:
        Xevs<-X.pca$sdev^2
        for(k in 1:ncol(X))
            if (sum(Xevs[1:k])/sum(Xevs) >= cutoff) break;
        Xpcs<-X %*% X.pca$rot[,1:k]
    }

    Xpcs
}

#calculate a permuttaion p-value matrix based on a permuted stat matrix:
PermPvs<-function(T0s){
    B=nrow(T0s); n=ncol(T0s)
    P0s<-matrix(1, nrow=B, ncol=n)
    for(j in 1:n)
        for(b in 1:B)
            P0s[b,j] = sum( abs(T0s[b,j]) < abs(T0s[-b,j]) )/(B-1)
    return(P0s)
}

## chkRank.drop.cols drops highly correlated variable for numerical problems
## return X if X has full column rank
chkRank.drop.cols <- function(X){

    p <- ncol(X)
    qr.X <- qr(X, tol = 1e-7, LAPACK = FALSE)
    rnkX <- qr.X$rank
    if (rnkX != p) {
       keep <- qr.X$pivot[seq_len(rnkX)]
       X <- X[, keep, drop = FALSE]
    }
  return(X = X)
}




# goal: transform data to repeated form for fitting function gee()
# dat: n by p data where sample size n and # of varibles
# k: number of traits
longForm <-function(dat, k){

  p <- ifelse(is.vector(dat), 1, ncol(dat))
  n <- ifelse(is.vector(dat), length(dat), nrow(dat))
  newdat <- matrix(0, ncol = p*k, nrow = n*k)
  for (r in 1:k) newdat[seq(r, n*k, k),(p*(r-1)+1):(p*r)] <- as.matrix(dat)
  return(newdat)
}

GEE <- function(traits, geno, Z = NULL,family = c("binomial", "gaussian"), corstr = "independence")
{
	 if (is.null(dim(traits))){
		 n <- length(traits);
		 k <- 1 } else {
    n <- nrow(traits)
	 k <- ncol(traits)}

    id <-rep(1:n, each = k)

    longX <- longForm(geno, k);
    longy <- matrix(as.vector( t(traits)), n*k, 1)
	 if (is.null(Z)) {longZ <- longForm(rep(1,n), k); mf <- formula(longy ~ 0+longZ)
		} else 	{longZ <- longForm(Z, k); mf <- formula(longy ~ longZ)}
  	 longZ <- chkRank.drop.cols(longZ)
	 fit0 <- gee(mf, id =id, family=family, corstr=corstr)
    R <- fit0$working.correlation
    muy <- fitted(fit0)
    res <- fit0$residuals

	 x <- cbind(longZ, longX)
	 Zind <- 1:ncol(longZ)
	 Xind <- (ncol(longZ)+1):ncol(x)

    U1 <- matrix(0, ncol=1, nrow=ncol(x))
    I <- matrix(0, ncol=ncol(x), nrow=ncol(x))

   ## First get the whole score vector and information matrix under null
    yy <- matrix(res, ncol = k, byrow = T)
    vy <- cov(yy,use="pairwise.complete.obs")
    U1<-matrix(0, ncol = 1,nrow = ncol(x))
    I <- matrix(0,ncol = ncol(x), nrow = ncol(x))

    if (family == "binomial")	v0 <- muy * (1 - muy)
	 if (family == "gaussian")	v0 <- rep(1, n*k)

    for (i in 1:n){

          muy1 <- muy[(k*(i-1)+1):(k*i)]
          v <- v0[(k*(i-1)+1):(k*i)]
   	    mux1 <-as.matrix(x[(k*(i-1)+1):(k*i),])
  	       dev  <- mux1*v
   	    yy <- res[(k*(i-1)+1):(k*i)]

		  # if k=num.traits is 1 or >1
          if (k > 1){
            V = diag(sqrt(c(v)))%*%R%*%diag(sqrt(c(v)))
            U0 = t(dev)%*%ginv(V)%*%yy
            I0 = t(dev)%*%ginv(V)%*%vy%*%t(ginv(V))%*%dev
          } else {
              V = v^2
              U0 = matrix(dev,ncol=1)*c(yy)/V
              I0 = matrix(dev,ncol=1)%*%matrix(dev,nrow=1)*c(vy)/V/V
          }


   	      U1 = U1+U0
   	      I = I+I0
   }

   ## calculate V where Us ~ N(0, V)
   ## under H0 the null distribution of the score vector for beta is asymptotically Normal
   ## I11 is the covariance matrix of gamma; I22 is the covariance matrix of beta;

	# recall followins defined
	# x <- cbind(longZ, longX)
	# Zind <- 1:ncol(longZ)
	# Xind <- (ncol(longZ)+1):ncol(x)

    Us = U1[Xind, ]
    I11 = as.matrix(I[Zind, Zind])
    I12 = as.matrix(I[Zind, Xind])
    I21 = as.matrix(I[Xind, Zind])
    I22 = as.matrix(I[Xind, Xind])
    dim(I21) = c(length(Xind), length(Zind))

    V <- I22 - (I21)%*%ginv(I11)%*%I12
    score <- t(Us)%*%ginv(V)%*% Us
    df <- ncol(longX)
    pval <- 1 - pchisq(score, df)
   return(list(U=Us,Cov=V, out=data.frame(df=df, Score=score,pval=pval)))
}





# Goal: obatin test statistics
Stat <- function(U, V, gamma)
{
	spu <- rep(NA,length(gamma))
	for (g in 1:length(gamma)) {
	  if (gamma[g]<Inf) spu[g]<-sum(U^gamma[g]) else spu[g]<-max(abs(U))
	}
	score <- sum((U %*% ginv(V)) *U)
	stat <- c(spu, score)
	return(stat)
}




# Goal: P values in SPU, aSPU, Score aSPU-Score
GEEspu.score <- function(U, V, gamma = c(1:8,Inf), B)
{

   spu.score <- Stat(U = U, V = V, gamma = gamma)

	p <- rep(NA, length(gamma)+1)
	T0s <- matrix(NA, nrow = B, ncol = length(gamma)+1)
   ### simulate the NULL
   U0 <- mvrnorm(B, rep(0,length(U)), V)
   for (b in 1:B){
	 T0s[b,] <- Stat(U = U0[b, ], V = V, gamma = gamma)
   }

   for (g in 1:ncol(T0s))  p[g] <- sum(abs(spu.score[g])<abs(T0s[,g]))/B
	P0s <- apply(T0s, 2, function(z) (B-rank(abs(z)))/(B-1))
	SPUminp0 <- apply(P0s[,1:length(gamma)],1,min)
	Paspu <- sum(SPUminp0 < min(p[1:length(gamma)]))/(B)

	SPUscoreminp0 <- apply(P0s,1,min)
	Paspu.score <- sum(SPUscoreminp0 < min(p))/(B)
	pvs <- c(head(p, length(gamma)), Paspu, tail(p, 1), Paspu.score)
   names(pvs)= c(paste("SPU", gamma, sep=""), "aSPU", "Score", "aSPU-Score")
	return(pvs)
}

MTaSPUsmallB <- function(Z, v, B, pow, transform = FALSE){
  # -- Z: matrix of summary Z-scores, SNPs in rows and traits in columns  
  # -- Or a vector of summary Z-scores for a single snp
  # -- v: output of estcov
  # --    estimated estimated correlation matrix  based on the summary Z-scores
  # -- tranform: if TRUE, the inference is made on transformed Z
  # -- B: number of Monte Carlo samples simulated to compute p-values 
  # -- results: compute p-values for SPU(gamma) i.e. pow=1:8, and infinity
  # --          aSPU, based on the minimum p-values over SPU(power)
  # --          each row for single SNP    
  
    
    if (transform){
        v <- ginv(v)
        Z <- tcrossprod(Z, v)
    }
    
    if (is.vector(Z)) {
        N <- 1; K <- length(Z)
    }	else {
        N <- dim(Z)[1]; K <- dim(Z)[2]
    }
 	 if (N ==1) Z <- matrix(Z, nrow=1)
	 dim(Z) = c(N, K)
    pval <- matrix(0, N, length(pow))
	 dim(pval) = c(N, length(pow))
    aSPU <- length(N)
    ponum <- pow[pow < Inf]
	 
    set.seed(1000)
    Z0 <- rmvnorm(B, mean = rep(0, nrow(v)), sigma = v)
	 ## SPU for integer power
    for(k0 in 1:length(ponum)){
        k <- ponum[k0]  
        z1 <- abs(rowSums(Z^k))
        z0b <- abs(rowSums(Z0^k))
        for(i in 1:N){
            pval[i,k0] <- (1+sum(z0b>z1[i]))/(B+1)
        }
    }
  
    ## SPU(max)
    if (Inf %in% pow){
        z1 <- rowMaxs(abs(Z))
        z0 <- rowMaxs(abs(Z0))
        for(i in 1:N){
			  pval[i,length(pow)] = (1+sum(z0>z1[i]))/(B+1)
		  }
    }
    
    ## aSPU
    p1m <- rowMins(pval)
    p0 <- matrix(NA, B, length(pow))
    for(k0 in 1:length(ponum)){
        k <- ponum[k0]
        zb <- abs(rowSums(Z0^k))
        p0[,k0] <- (1+B-rank(abs(zb)))/B
    }

    
    if (Inf %in% pow){
      zb <- rowMaxs(abs(Z0))
      p0[,length(pow)] = (1+B-rank(zb))/B
    }
    p0m = rowMins(p0)
    for(i in 1:N){
        aSPU[i] = (1+sum(p0m<p1m[i]))/(B+1)
    }
	 pval <- cbind(pval, aSPU)
    if(Inf %in% pow) s <- c(paste("SPU", ponum, sep=""),"SPUInf","MTaSPUs") else {
      s <- c(paste("SPU", ponum, sep=""),"MTaSPUs")}
    colnames(pval) <- s
    rownames(pval) <- rownames(Z)
    return(pval)
}











MTaSPUsB1e8 <- function(Z, v, pow, transform = FALSE){
   # -- B: number of Monte Carlo simulation is fixed at B=1e8
   # -- Z: matrix of summary Z-scores, SNPs in rows and traits in columns  
   # -- Or a vector of summary Z-scores for a single snp
   # -- v: output of estcov
   # --    estimated estimated correlation matrix  based on the summary Z-scores
   # -- tranform: if TRUE, the inference is made on transformed Z
   # -- results: compute p-values for SPU(gamma) i.e. pow=1:8, and infinity
   # --          aSPU, based on the minimum p-values over SPU(power)
   # --          each row for single SNP    
  
    
   B<- 1e7; B2<- 10
	if (transform){
	    v <- ginv(v)
	    Z <- tcrossprod(Z, v)
	}
    
	if (is.vector(Z)) {
		N <- 1; K <- length(Z)
	} else {
		N <- dim(Z)[1]; K <- dim(Z)[2]
	}
	if (N ==1) Z <- matrix(Z, nrow=1)
	ponum <- pow[pow < Inf]
	Zpu = apply(Z, 1, function(zj){
	  tmp = rep(NA,length(pow))
	  for(k0 in 1:length(ponum)){
	    tmp[k0]<- abs(sum(zj^ponum[k0]))
	  }
	  if (Inf %in% pow){tmp[length(pow)] = max(abs(zj))}
	  return(tmp)
	})

	P1perm <- matrix(0, N, length(pow))
	Qq <- matrix(0, 50, length(pow))
	for(b in 1:B2){
	  zb <- rmvnorm(B, sigma = v)
	  for(k0 in 1:length(ponum)){
		 k <- ponum[k0]  
	    res = abs(rowSums(zb^k))
	    for(j in 1:N){
	      P1perm[j,k0] = P1perm[j,k0]+mean(res>Zpu[k0,j])/B2
	    }
	    res = c(res, Qq[,k0])
	    Qq[,k0] = -sort.int(-res,partial=50)[1:50]
	  }
	  if (Inf %in% pow){
		 res = rowMaxs(abs(zb))
		 for(i in 1:N){
		 	P1perm[i,length(pow)] = P1perm[i,length(pow)]+mean(res>Zpu[length(pow),i])/B2
		 }
		 res = c(res, Qq[,length(pow)])
		 Qq[,length(pow)] = -sort.int(-res,partial=50)[1:50]
	  }
	}


	## aSPU
	Qq = apply(Qq, 2, sort)
	Q0 = rbind(0,Qq)
	pr0 = c(1,50:1/(B*B2))

	pvi = matrix(NA, B,length(pow))
	Pmin = rowMins(P1perm)
	aSPU = numeric(N)
	for(b in 1:B2){
	  zb = rmvnorm(B, sigma = v)
	  for(k0 in 1:length(ponum)){
		 k <- ponum[k0]  
	    pvi[,k0] = approx(Q0[,k0], pr0, abs(rowSums(zb^k)), method='linear',rule=2)$y
	  }
	  if (Inf %in% pow){
	  	pvi[,length(pow)] = approx(Q0[,length(pow)], pr0, rowMaxs(abs(zb)), method='linear',rule=2)$y
	  }
  	  tmp = rowMins(pvi)
	  for(i in 1:N){
	    aSPU[i] = aSPU[i] + mean(tmp < Pmin[i])/B2
	  }
	}
	P1perm = cbind(P1perm, aSPU)
	P1perm[P1perm==0] = 1/(B*B2)
   if(Inf %in% pow) s <- c(paste("SPU", ponum, sep=""),"SPUInf","MTaSPUs") else {
      s <- c(paste("SPU", ponum, sep=""),"MTaSPUs")}
	colnames(P1perm) <- s
   rownames(P1perm) <- rownames(Z)
	return(P1perm)
}
