pbad2way<-function(formula, data, est = "mom", nboot = 599){
  
  if (missing(data)) {
    mf <- model.frame(formula)
  } else {
    mf <- model.frame(formula, data)
  }
  cl <- match.call()
  est <- match.arg(est, c("mom", "onestep", "median"), several.ok = FALSE)
  
  J <- nlevels(mf[,2])
  K <- nlevels(mf[,3])
  
  alpha=.05
  conall = TRUE 
  op=FALSE
  pro.dis=FALSE
  MM=FALSE
  
  grp <- NA
  JK <- J * K
  est <- get(est)
  
  nfac <- tapply(mf[,1], list(mf[,2],mf[,3]), length, simplify = FALSE)
  nfac1 <- nfac[unique(mf[,2]), unique(mf[,3])]    ## reordering factor levels
  #data$row <- unlist(alply(nfac1, 1, sequence), use.names = FALSE)
  mf$row <-  unlist(alply(nfac1, 1, sequence), use.names = FALSE)
  
  dataMelt <- melt(mf, id = c("row", colnames(mf)[2], colnames(mf)[3]), measured = mf[,1])
  dataWide <- cast(dataMelt, as.formula(paste(colnames(dataMelt)[1], "~", colnames(mf)[2], "+", colnames(mf)[3]))) 
  dataWide$row <- NULL
  #dataMelt <- melt(data, id = c("row", colnames(mf)[2], colnames(mf)[3]), measured = mf[,1])
  #dataWide <- cast(dataMelt, as.formula(paste(colnames(dataMelt)[1], "~", colnames(mf)[2], "+", colnames(mf)[3]))) 
  #dataWide$row <- NULL
  x <- dataWide
  
  if(is.matrix(x)) x <- as.data.frame(x)
  x <- listm(x)
  #if(!is.na(grp[1])) {
  #  yy <- x
  #  for(j in 1:length(grp))
  #    x[[j]] <- yy[[grp[j]]]
  #}
  #if(!is.list(x)) stop("Data must be stored in list mode or a matrix.")
  for(j in 1:JK) {
    xx <- x[[j]]
    x[[j]] <- xx[!is.na(xx)]
  }
  #
  # Create the three contrast matrices
  #
  if(!conall){
    ij <- matrix(c(rep(1, J)), 1, J)
    ik <- matrix(c(rep(1, K)), 1, K)
    jm1 <- J - 1
    cj <- diag(1, jm1, J)
    for(i in 1:jm1)
      cj[i, i + 1] <- 0 - 1
    km1 <- K - 1
    ck <- diag(1, km1, K)
    for(i in 1:km1)
      ck[i, i + 1] <- 0 - 1
    conA <- t(kron(cj, ik))
    conB <- t(kron(ij, ck))
    conAB <- t(kron(cj, ck))
    conAB <- t(kron(abs(cj), ck))
  }
  if(conall){
    temp<-con2way(J,K)
    conA<-temp$conA
    conB<-temp$conB
    conAB<-temp$conAB
  }
  ncon <- max(nrow(conA), nrow(conB), nrow(conAB))
  if(!is.na(grp[1])){  # Only analyze specified groups.
    xx<-list()
    for(i in 1:length(grp))xx[[i]]<-x[[grp[i]]]
    x<-xx
  }
  mvec<-NA
  for(j in 1:JK){
    temp<-x[[j]]
    temp<-temp[!is.na(temp)] # Remove missing values.
    x[[j]]<-temp
    mvec[j]<-est(temp)
  }
  bvec<-matrix(NA,nrow=JK,ncol=nboot)
  #set.seed(2) # set seed of random number generator so that
  #             results can be duplicated.
  #print("Taking bootstrap samples. Please wait.")
  for(j in 1:JK){
    data<-matrix(sample(x[[j]],size=length(x[[j]])*nboot,replace=TRUE),nrow=nboot)
    bvec[j,]<-apply(data,1,est) # J by nboot matrix, jth row contains bootstrapped  estimates for jth group
    naind <- which(is.na(bvec[j,]))
    if (length(naind) > 0) bvec[j,][naind] <- mean(bvec[j,], na.rm = TRUE)    ## fix for missing est values
  }
  
  bconA<-t(conA)%*%bvec #C by nboot matrix
  tvecA<-t(conA)%*%mvec
  tvecA<-tvecA[,1]
  tempcenA<-apply(bconA,1,mean)
  veczA<-rep(0,ncol(conA))
  bconA<-t(bconA)
  smatA<-var(bconA-tempcenA+tvecA)
  bconA<-rbind(bconA,veczA)
  if(!pro.dis){
    if(!op) dv<-mahalanobis(bconA,tvecA,smatA)
    if(op){
      dv<-out(bconA)$dis
    }}
  if(pro.dis)dv=pdis(bconA,MM=MM)
  bplus<-nboot+1
  sig.levelA<-1-sum(dv[bplus]>=dv[1:nboot])/nboot
  bconB<-t(conB)%*%bvec #C by nboot matrix
  tvecB<-t(conB)%*%mvec
  tvecB<-tvecB[,1]
  tempcenB<-apply(bconB,1,mean)
  veczB<-rep(0,ncol(conB))
  bconB<-t(bconB)
  smatB<-var(bconB-tempcenB+tvecB)
  bconB<-rbind(bconB,veczB)
  if(!pro.dis){
    if(!op)dv<-mahalanobis(bconB,tvecB,smatB)
    if(op){
      dv<-out(bconA)$dis
    }}
  if(pro.dis)dv=pdis(bconB,MM=MM)
  sig.levelB<-1-sum(dv[bplus]>=dv[1:nboot])/nboot
  bconAB<-t(conAB)%*%bvec #C by nboot matrix
  tvecAB<-t(conAB)%*%mvec
  tvecAB<-tvecAB[,1]
  tempcenAB<-apply(bconAB,1,mean)
  veczAB<-rep(0,ncol(conAB))
  bconAB<-t(bconAB)
  smatAB<-var(bconAB-tempcenAB+tvecAB)
  bconAB<-rbind(bconAB,veczAB)
  if(!pro.dis){
    if(!op)dv<-mahalanobis(bconAB,tvecAB,smatAB)
    if(op){
      dv<-out(bconAB)$dis
    }}
  if(pro.dis)dv=pdis(bconAB,MM=MM)
  sig.levelAB<-1-sum(dv[bplus]>=dv[1:nboot])/nboot
  
  result <- list(Qa = NA, A.p.value=sig.levelA, Qb=NA, B.p.value=sig.levelB, Qab=NA, AB.p.value=sig.levelAB, call = cl, varnames = colnames(mf))
  class(result) <- c("t2way")
  result
}
