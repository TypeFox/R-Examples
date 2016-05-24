# nearest neighbor pooled linear dimension reduction according to
# Hastie and Tibshirani,
# IEEE Trans. Pattern Analysis and Machine Intelligence 18 (1996), 607-616
ncoord <- function(xd, clvecd, nn=50, weighted=FALSE,
                    sphere="mcd", orderall=TRUE, countmode=1000, ...){
  z <- x <- as.matrix(xd)
  if (is.matrix(sphere)){
    cv <- sphere
    sphere <- "matrix"
  }
  if (sphere!="none"){
#  require(MASS)
    if (sphere=="matrix")
      Sig <- cv
    else
      Sig <- cov.rob(x, method=sphere, nsamp=500)$cov
    Tds <- tdecomp(Sig)
    Y <- solve(Tds)
    z <- x %*% Y
  }
  clvec <- as.integer(clvecd)
  clf <- factor(clvec)
  cll <- as.integer(levels(clf)) 
  cln <- length(cll)
  p <- ncol(z)
  n <- nrow(z)
  B <- matrix(0,ncol=p,nrow=p)
  for (i in 1:n){
    if(countmode*round(i/countmode)==i)
      cat("Processing point ",i," of ",n,"\n")
    Bi <- matrix(0,ncol=p,nrow=p)
    za <- sweep(z,2,z[i,])
    mds <- rowSums(za*za)
    if (orderall)
      omah <- order(mds)[1:nn]
    else{
      omah <- c()
      maxmds <- max(mds)
      for (j in 1:nn){
        argmin <- which.min(mds)
        omah <- c(omah,argmin)
        mds[argmin] <- maxmds+1
      }
    }
    mi <- colMeans(z[omah,])
    for (j in 1:cln){
      nj <- sum(clvec[omah]==cll[j])
      pij <- nj/nn
      v <- if (pij==0) rep(0,p)
           else colMeans(z[omah,][clvec[omah]==cll[j],,drop=FALSE])-mi
      Bi <- Bi+ pij*(v %*% rbind(v))
    }
    if (weighted){
      sb <- sum(diag(Bi))
      if (sb>0)
        Bi <- Bi/sb
    }
    B <- B+Bi
  }
  B <- B/n
  em <- eigen(B, symmetric=TRUE)
  units <- em$vectors
  if (sphere!="none")
    units <- Y %*% units
  proj <- x %*% units 
  list(ev=em$values, units=units, proj=proj)
}

# asymmetric robust
# nearest neighbor pooled linear dimension reduction 
ancoord <- function(xd, clvecd, clnum=1, nn=50, method="mcd",
                    countmode=1000, ...){
#  require(MASS)
  x <- as.matrix(xd)
  p <- ncol(x)
  n <- nrow(x)
  clvec <- as.integer(clvecd)
  dcl <- as.integer(clnum)
  ci <- clvec==dcl
  clxf <- x[ci,]
  nc <- sum(ci)
  quant <- min(floor(3*(nrow(clxf) + ncol(clxf) + 1)/4),nrow(clxf)-2)
  cv <- cov.rob(clxf,quantile.used=quant,method=method,nsamp=500)
  S1 <- cv$cov
  cinv <- solvecov(S1)$inv
  B <- matrix(0,ncol=p,nrow=p)
  w <- 0
  repeat{
    for (i in 1:nc){
      if(countmode*round(i/countmode)==i)
        cat("Processing point ",i," of ",nc,"\n")
      Bi <- matrix(0,ncol=p,nrow=p)
      mds <- mahalanobis(x,center=clxf[i,],cov=cinv,inverted=TRUE)
      omah <- order(mds)[1:nn]
      wi <- 1    
      mi <- colMeans(x[omah,])
      ni <- sum(clvec[omah]==dcl)
      nr <- nn-ni
      wi <- ni*nr
      if (wi>0){
        vi <- colMeans(x[omah,][ci[omah],,drop=FALSE])-mi
        vr <- colMeans(x[omah,][!ci[omah],,drop=FALSE])-mi
        Bi <- Bi+ ni*(vi %*% rbind(vi))+ nr*(vr %*% rbind(vr))
      }
      sb <- sum(diag(Bi))
      if (sb>0)
        Bi <- Bi/(nn*sb)
      B <- B+Bi
    }
    if (!identical(B,matrix(0,ncol=p,nrow=p)))
      break
    else{
      if (nn<nc+1)
        nn <- nc+1
      else{
        warning("Estimated between groups matrix is zero!")
        break
      }
    }
  }
  Tm <- tdecomp(S1)
  Tinv <- solve(Tm)
  Z <- t(Tinv) %*% B %*% Tinv
  dc <- eigen(Z, symmetric=TRUE)
  units <- Tinv %*% dc$vectors
  proj <- x %*% units    
  list(ev=dc$values, units=units, proj=proj, nn=nn)
}

# quadratic dimension reduction according to Young, Marco and Odell,
# Journal Stat. Plann. Inf. 17 (1986), 307-319; computation according to
# Roehl and Weihs, in Gaul & Locarek-Junge (1999), 253.
mvdcoord <- function(xd, clvecd, clnum=1, sphere="mcd", ...){
  x <- as.matrix(xd)
  if (is.matrix(sphere)){
    cv <- sphere
    sphere="matrix"
  }
#  require(MASS)
  if (sphere!="none"){
    if (sphere=="matrix")
      Sig <- cv
    else
      Sig <- cov.rob(x, method=sphere, nsamp=500)$cov
    Tds <- tdecomp(Sig)
    Y <- solve(Tds)
    z <- x %*% Y
  }
  else
    z <- x
  clvec <- as.integer(clvecd)
  clf <- factor(clvec)
  cll <- as.integer(levels(clf)) 
  clnum <- length(cll)
  p <- ncol(z)
  mx <- vx <- list()
  for (i in 1:clnum){
    mx[[i]] <- colMeans(z[clvecd==cll[i],])
    vx[[i]] <- cov(z[clvecd==cll[i],])
  }
  meandiff <- vardiff <- c()
  for (i in 2:clnum){
    meandiff <- cbind(meandiff,mx[[i]]-mx[[1]])
    vardiff <- cbind(vardiff,vx[[i]]-vx[[1]])
  }
  M <- cbind(meandiff,vardiff)
  em <- eigen(M %*% t(M), symmetric=TRUE)
  units <- em$vectors
  if (sphere!="none")
    units <- Y %*% units
  proj <- x %*% units 
  list(ev=em$values, units=units, proj=proj)
}

# mahalanodisc=vector of mahalanobis distances from n1 points of x1
# to n-n1 points from x2
# modus: see mahal in robcoord
mahalanodisc <- function (x2, mg, covg, modus="square") {
  covinv <- solvecov(covg)$inv
  md <- switch(modus,
         md=sqrt(mahalanobis(x2,mg,covinv,inverted=TRUE)),
         mahalanobis(x2,mg,covinv,inverted=TRUE))
  md
}
# dist:n-n1 Mahalanobis distances, mg: mean(x1), covg: Covariance(x1)

# weight function for robcoord
cweight <- function(x,ca){
  out <- 1
  if (x > ca)
    out <- ca/x
  out
}

adcoord <- function(xd, clvecd, clnum=1) {
  x <- as.matrix(xd)
  clvec <- as.integer(clvecd)
  dcl <- as.integer(clnum)
  ci <- clvec==dcl
  n <- nrow(x)
  p <- ncol(x)
  cln <- sum(ci)
  clx <- rep(0, times=p*cln)
  clxc <- rep(0, times=p*(n-cln))
  dim(clx) <- c(cln,p)
  dim(clxc) <- c((n-cln),p)
  for (j in 1:p){
    clx[,j] <- x[,j][ci]
    clxc[,j] <- x[,j][!ci]
  }
  S1 <- cov(clx)
  S2 <- cov(clxc)
  W1 <- (cln-1)*S1
  S <- cov(x)
  B <- (n*(n-1)*S - cln*W1 - (n-cln)*(n-cln-1)*S2)/(n-cln)
  Tm <- tdecomp(S1)
  Tinv <- solve(Tm)
  Z <- t(Tinv) %*% B %*% Tinv
  dc <- eigen(Z, symmetric=TRUE)
  units <- Tinv %*% dc$vectors
  proj <- x %*% units    
  list(ev=dc$values, units=units, proj=proj)
}

# "robustifizierte" 1-Cluster-Diskriminanzkoordinaten (durchschnittlicher
# Innerhalb-Abstand vs. Abstand nach ausserhalb, letzterer gewichtet
# mit c(Mahal gross)/Mahal(x_j-mean1).
# Projektionen, Eigenwerte
# x: Daten, clvec: Clusterindikatorvektor, clnum: Nummer des
# zu trennenden Clusters ,
# mahal="square": Squared Mahalanobis distance is used
# mahal="md": Mahalanobis distance is used
# subsample: size of subsample of cluster to use (0=all)
# countmode=output of current point number
awcoord <- function(xd, clvecd, clnum=1, mahal="square", method="classical",
                     clweight=switch(method,classical=FALSE,TRUE), alpha=0.99,
                     subsample=0, countmode=1000, ...) {
  x <- as.matrix(xd)
#  require(MASS)
  n <- nrow(x)
  p <- ncol(x)
  dcl <- as.integer(clnum)
  clfull <- as.integer(clvecd)
  cln <- sum(clfull==dcl)
  if (subsample==0){
    clvec <- as.integer(clfull==dcl)
    cn <- cln
    subs <- NULL
    clxf <- clx <- x[clvec==1,]
  }
  else{
    subs <- sample((1:n)[clfull==dcl],subsample)
    clvec <- 2*(clfull==dcl)
    clvec[subs] <- 1
    clxf <- x[clvec==2,]
    clx <- x[clvec==1,]
    cn <- subsample
  }
  clxc <- x[clvec==0,]
  clxa <- if (clweight) x else clxc
  quant <- min(floor(3*(nrow(clxf) + ncol(clxf) + 1)/4),nrow(clxf)-2)
  cv <- cov.rob(clxf, quantile.used=quant,method=method,nsamp=500)
  S1 <- cv$cov
  mg <- cv$center
  mah <- mahalanodisc(clxa, mg, S1, modus=mahal)
  wg <- switch(mahal,
     md=sapply(mah,cweight,sqrt(qchisq(alpha,p))),
     sapply(mah,cweight,qchisq(alpha,p))) 
   if (clweight){
     wg0 <- wg[clvec==0]
     wg1 <- wg[clvec==1]
     wsum0 <- sum(wg0)
     wsum1 <- sum(wg1)
   }
  else{
    wg0 <- wg
    wsum0 <- sum(wg)
  }
   d1 <- d2 <- rep(0,p)
   D12 <- D22 <- matrix(0,ncol=p, nrow=p)
   for(i in 1:cn){
     if (clweight)
       D12 <- D12 + wg1[i]*clx[i,] %*% t(clx[i,])
     else
       D12 <- D12 + clx[i,] %*% t(clx[i,])
     if (clweight)
       d1 <- d1+wg1[i]*clx[i,]
     else
       d1 <- d1+clx[i,]
   }
   for (j in 1:(n-cln)){
     D22 <- D22 + wg0[j]*clxc[j,] %*% t(clxc[j,])
     d2 <- d2+wg0[j]*clxc[j,]
   }
   D21 <- d1 %*% t(d2)
   if (clweight)
     B <- wsum0*D12-D21-t(D21)+wsum1*D22
   else
     B <- wsum0*D12-D21-t(D21)+cn*D22
  if (clweight)
    wsum <- sum(outer(wg0,wg1))
  else              
     wsum <- wsum0*cn
  B <- B/wsum
  Tm <- tdecomp(S1)
  Tinv <- solve(Tm)
  Z <- t(Tinv) %*% B %*% Tinv
  dc <- eigen(Z, symmetric=TRUE)
  units <- Tinv %*% dc$vectors
  proj <- x %*% units    
  list(ev=dc$values, units=units, proj=proj, wg=wg)
}

