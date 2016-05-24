 library(Matrix) # fuer rankMatrix
 library(MASS) # fuer ginv

 allelesperlocus <- function(vec){
   Nralleles <- sort(unique(vec))
   alleles=Nralleles
 }

 nrallelesperlocus <- function(vec){
   Nralleles <- sort(unique(vec))
   nralleles=length(Nralleles)
 }

 create.WandB <- function(vec){
   elements <- sort(unique(vec))
   res <- t(sapply(vec,function(y){y==elements}))
   Wk <- matrix(0,dim(res)[1],dim(res)[2])
   Wk[res] <- 1
   Bk <- matrix(1,dim(res)[1],dim(res)[2])
   Bk[is.na(res)] <- 0
   list(Wk,   Bk)
 }

 createM.step2 <- function(locimat,forQ,forNralleles,formini,formaxi,formini_m,formaxi_m,forV,forW,forB,foruse_m1,forpart2a,forpart2b){
   r <- forNralleles[locimat[1]]
   ncols <- ncol(forQ)
   nrows <- r-1
   H1_2 <- matrix(0,nrows,ncols)
   H2_2 <- matrix(0,nrows,nrows)
   diag(H2_2) <- 1
   H3_2 <- rep(-1,nrows)
   Hq_2 <- cbind(H1_2,H2_2,H3_2)
   tHq_2 <- t(Hq_2)
   use_m2 <- c(1:nrow(forB))[forB[,formini[locimat[1]]]==1]
   use <- intersect(foruse_m1,use_m2)
   Vq_2 <- (forV[foruse_m1,use_m2])
   ginvVq_3 <- solve(forV[use_m2,use_m2])
   Xq_2 <- as.matrix(cbind(forQ,forW[,formini[locimat[1]]:formaxi[locimat[1]]]))
   Xqred_2 <- Xq_2[forB[,formini[locimat[1]]]==1,]
   part4 <- ginv(crossprod(Xqred_2,ginvVq_3)%*%Xqred_2)
   part5 <- forpart2a%*%(Vq_2%*%ginvVq_3%*%Xqred_2)
   partM <- unlist(forpart2b%*%part5%*%tcrossprod(part4,Hq_2))
   partM
 }

 calculate.threshold <- function(alleles,forp0,forM){
   sel <- (alleles[2]):(alleles[3])
   p1 <- forp0[sel]
   p2 <- forp0[sel,]
   p3 <- solve(forM[sel,sel])
   T <- crossprod(p1,p3)%*%p2
   1-pchisq(T,alleles[1]-1)
 }

 gwerAM <- function(A,Q,W,varianceA,varianceE,Nrsim=10000,k=c(0,1,2,5),alpha=c(0.05)) {
 
   Alleles <- apply(W,2,allelesperlocus)
   Nralleles <- apply(W,2,nrallelesperlocus)
   Nrgenotypes <- nrow(W)
   Nrloci <- length(Nralleles)

   out <- apply(W,2,create.WandB)
   Wstart <- lapply(out,"[[",1)
   Bstart <- lapply(out,"[[",2)
   W2 <- matrix(unlist(Wstart),nrow=Nrgenotypes) 
   B <- matrix(unlist(Bstart),nrow=Nrgenotypes) 
   Q <- cbind(rep(1,Nrgenotypes),Q)

   Rmatrix <- matrix(0,Nrgenotypes,Nrgenotypes)
   diag(Rmatrix) <-1
   Partr <- Rmatrix*varianceE
   Parta <- A*varianceA

   V <- as.matrix(Partr+Parta)
   maxi <- cumsum(Nralleles)
   mini <- c(1,(maxi+1)[1:(Nrloci-1)])

   M <- matrix(NA,sum(Nralleles-1),sum(Nralleles-1))
   maxi_m <- maxi-c(1:Nrloci)
   mini_m <- mini+(-c(1:Nrloci)+1)
 
   for(i in 1:Nrloci){
     r <- Nralleles[i]
     ncols <- ncol(Q)
     nrows <- r-1
     H1_1 <- matrix(0,nrows,ncols)
     H2_1 <- matrix(0,nrows,nrows)
     diag(H2_1) <- 1
     H3_1 <- rep(-1,nrows)
     Hq_1 <- cbind(H1_1,H2_1,H3_1)
     tHq_1 <- t(Hq_1)
     Vq_1 <- (V[B[,mini[i]]==1,B[,mini[i]]==1])
     Invtest_1 <- solve(Vq_1)
     Xq_1 <- as.matrix(cbind(Q,W2[,mini[i]:maxi[i]]))
     Xqred_1 <- Xq_1[B[,mini[i]]==1,]
     part2a <- crossprod(Xqred_1,Invtest_1)
     part2 <- ginv(part2a%*%Xqred_1)
     part2b <- Hq_1%*%part2
     use_m1 <- c(1:nrow(B))[B[,mini[i]]==1]
     partMsame <- Hq_1%*%tcrossprod(part2,Hq_1)
     M[(mini_m[i]):(maxi_m[i]),(mini_m[i]):(maxi_m[i])] <- partMsame
     if(i!=Nrloci){
       locus2 <- cbind(c((i+1):Nrloci),c((i+1):Nrloci))
       partMdiff <- apply(
           locus2,
           1,
           createM.step2,
           forQ=Q,
           forNralleles=Nralleles,
           formini=mini,
           formaxi=maxi,
           formini_m=mini_m,
           formaxi_m=maxi_m,
           forV=V,
           forW=W2,
           forB=B,
           foruse_m1=use_m1,
           forpart2a=part2a,
           forpart2b=part2b
           )
       keep <- matrix(unlist(partMdiff),nrow=nrow(partMsame))
       M[(mini_m[i]):(maxi_m[i]),(mini_m[i+1]):(maxi_m[Nrloci])] <- keep
       M[(mini_m[i+1]):(maxi_m[Nrloci]),(mini_m[i]):(maxi_m[i])] <- t(keep)
     }
   }#i

   mysvd <- svd(M)
   le <- mysvd$d
   u <- mysvd$u
   P <-  u%*%diag(sqrt(le))
   psim <- matrix(NA,Nrsim,length(k))
   rankP <- rankMatrix(P)
   innersim <- matrix(rnorm(Nrsim*rankP,mean=0,sd=1),rankP,Nrsim)

   step2 <- cbind(Nralleles,mini_m,maxi_m)
   for(i in 1:Nrsim){
     p0 <- P%*%innersim[,i]
     pvalue <- apply(step2,1,calculate.threshold,forp0=p0,forM=M)
     psim[i,] <- sort(pvalue)[k+1]
   }
   colnames(psim) <- paste("k=",k,sep="")
   apply(psim,2,quantile,probs=alpha,na.rm=TRUE)
 }

