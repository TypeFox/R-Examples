sppba <- function(formula, id, data, est = "mom", avg = TRUE, nboot = 500, MDIS = FALSE){
  
  if (missing(data)) {
    mf <- model.frame(formula)
  } else {
    mf <- model.frame(formula, data)
  }
  cl <- match.call()
  est <- match.arg(est, c("mom", "onestep", "median"), several.ok = FALSE)
  
 
  mf1 <- match.call()
  m <- match(c("formula", "data", "id"), names(mf1), 0L)
  mf1 <- mf1[c(1L, m)]
  mf1$drop.unused.levels <- TRUE
  mf1[[1L]] <- quote(stats::model.frame)
  mf1 <- eval(mf1, parent.frame())  

  random1 <- mf1[, "(id)"]
  depvar <- colnames(mf)[1]

## check which one is the within subjects factor
  if (all(length(table(random1)) == table(mf[,3]))) {
    ranvar <- colnames(mf)[3]
    fixvar <- colnames(mf)[2]
  } else {
    ranvar <- colnames(mf)[2]
    fixvar <- colnames(mf)[3]
  }

  MC <- FALSE
  K <- length(table(mf[, ranvar]))  ## number of repeated measurements
  J <- length(table(mf[, fixvar]))  ## number of levels
  p <- J*K
  grp <- 1:p
  est <- get(est)  

  fixsplit <- split(mf[,depvar], mf[,fixvar])
  indsplit <- split(mf[,ranvar], mf[,fixvar])
  dattemp <- mapply(split, fixsplit, indsplit, SIMPLIFY = FALSE)
  data <- do.call(c, dattemp)
  x <- data

  jp<-1-K
  kv<-0
  kv2<-0
  for(j in 1:J){
    jp<-jp+K
    xmat<-matrix(NA,ncol=K,nrow=length(x[[jp]]))
    for(k in 1:K){
      kv<-kv+1
      xmat[,k]<-x[[kv]]
    }
    xmat<-elimna(xmat)
    for(k in 1:K){
      kv2<-kv2+1
      x[[kv2]]<-xmat[,k]
    }
  }
  xx<-x
  #if(SEED)set.seed(2) # set seed of random number generator so that
  #             results can be duplicated.
  # Next determine the n_j values
  nvec<-NA
  jp<-1-K
  for(j in 1:J){
    jp<-jp+K
    nvec[j]<-length(x[[jp]])
  }
  #
  # Now take bootstrap samples from jth level
  # of Factor A.
  #
  bloc<-matrix(NA,nrow=J,ncol=nboot)
  #print("Taking bootstrap samples. Please wait.")
  mvec<-NA
  ik<-0
  
  for(j in 1:J){
    #paste("Working on level ",j," of Factor A")
    x<-matrix(NA,nrow=nvec[j],ncol=K)
    #
    for(k in 1:K){
      ik<-ik+1
      x[,k]<-xx[[ik]]
      if(!avg)mvec[ik]<-est(xx[[ik]])
    }
    tempv<-apply(x,2,est)
    data<-matrix(sample(nvec[j],size=nvec[j]*nboot,replace=TRUE),nrow=nboot)
    bvec<-matrix(NA,ncol=K,nrow=nboot)
    for(k in 1:K){
      temp<-x[,k]
      bvec[,k]<-apply(data,1,rmanogsub,temp,est) # An nboot by K matrix
    }
    if(avg){
      mvec[j]<-mean(tempv)
      bloc[j,]<-apply(bvec,1,mean)
    }
    if(!avg){
      if(j==1)bloc<-bvec
      if(j>1)bloc<-cbind(bloc,bvec)
    }
  }
  
  
  if(avg){
    d<-(J^2-J)/2
    con<-matrix(0,J,d)
    id<-0
    Jm<-J-1
    for (j in 1:Jm){
      jp<-j+1
      for(k in jp:J){
        id<-id+1
        con[j,id]<-1
        con[k,id]<-0-1
      }}}
  if(!avg){
    MJK<-K*(J^2-J)/2 # NUMBER OF COMPARISONS
    JK<-J*K
    MJ<-(J^2-J)/2
    cont<-matrix(0,nrow=J,ncol=MJ)
    ic<-0
    for(j in 1:J){
      for(jj in 1:J){
        if(j<jj){
          ic<-ic+1
          cont[j,ic]<-1
          cont[jj,ic]<-0-1
        }}}
    tempv<-matrix(0,nrow=K-1,ncol=MJ)
    con1<-rbind(cont[1,],tempv)
    for(j in 2:J){
      con2<-rbind(cont[j,],tempv)
      con1<-rbind(con1,con2)
    }
    con<-con1
    if(K>1){
      for(k in 2:K){
        con1<-push(con1)
        con<-cbind(con,con1)
      }}}
  
  if(!avg)bcon<-t(con)%*%t(bloc) #C by nboot matrix
  if(avg)bcon<-t(con)%*%(bloc)
  tvec<-t(con)%*%mvec
  tvec<-tvec[,1]
  tempcen<-apply(bcon,1,mean)
  vecz<-rep(0,ncol(con))
  bcon<-t(bcon)
  temp=bcon
  for(ib in 1:nrow(temp))temp[ib,]=temp[ib,]-tempcen+tvec
  bcon<-rbind(bcon,vecz)
  if(!MDIS){
    if(!MC)dv=pdis(bcon,center=tvec)
    #if(MC)dv=pdisMC(bcon,center=tvec)
  }
  if(MDIS){
    smat<-var(temp)
    bcon<-rbind(bcon,vecz)
    chkrank<-qr(smat)$rank
    if(chkrank==ncol(smat))dv<-mahalanobis(bcon,tvec,smat)
    if(chkrank<ncol(smat)){
      smat<-ginv(smat)
      dv<-mahalanobis(bcon,tvec,smat,inverted=T)
    }}
  bplus<-nboot+1
  sig.level<-1-sum(dv[bplus]>=dv[1:nboot])/nboot
  
  result <- list(test = tvec, p.value = sig.level, call = cl)
  class(result) <- c("spp")
  result
}
