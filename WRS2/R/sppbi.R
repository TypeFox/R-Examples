sppbi <- function(formula, id, data, est = "mom", nboot = 500){
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
  
  JK<-J*K
  MJ<-(J^2-J)/2
  MK<-(K^2-K)/2
  JMK<-J*MK
  Jm<-J-1
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
    }}
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
  # of Factor A and average K  corresponding estimates
  # of location.
  #
  bloc<-matrix(NA,ncol=J,nrow=nboot)
  #print("Taking bootstrap samples. Please wait.")
  mvec<-NA
  it<-0
  for(j in 1:J){
    #paste("Working on level ",j," of Factor A")
    x<-matrix(NA,nrow=nvec[j],ncol=MK)
    #
    im<-0
    for(k in 1:K){
      for(kk in 1:K){
        if(k<kk){
          im<-im+1
          kp<-j*K+k-K
          kpp<-j*K+kk-K
          x[,im]<-xx[[kp]]-xx[[kpp]]
          it<-it+1
          mvec[it]<-est(x[,im])
        }}}
    data<-matrix(sample(nvec[j],size=nvec[j]*nboot,replace=TRUE),nrow=nboot)
    bvec<-matrix(NA,ncol=MK,nrow=nboot)
    mat<-listm(x)
    for(k in 1:MK){
      temp<-x[,k]
      bvec[,k]<-apply(data,1,rmanogsub,temp,est) # An nboot by MK matrix
    }
    if(j==1)bloc<-bvec
    if(j>1)bloc<-cbind(bloc,bvec)
  }
  #
  MJMK<-MJ*MK
  con<-matrix(0,nrow=JMK,ncol=MJMK)
  cont<-matrix(0,nrow=J,ncol=MJ)
  ic<-0
  for(j in 1:J){
    for(jj in 1:J){
      if(j<jj){
        ic<-ic+1
        cont[j,ic]<-1
        cont[jj,ic]<-0-1
      }}}
  tempv<-matrix(0,nrow=MK-1,ncol=MJ)
  con1<-rbind(cont[1,],tempv)
  for(j in 2:J){
    con2<-rbind(cont[j,],tempv)
    con1<-rbind(con1,con2)
  }
  con<-con1
  if(MK>1){
    for(k in 2:MK){
      con1<-push(con1)
      con<-cbind(con,con1)
    }}
  bcon<-t(con)%*%t(bloc) #C by nboot matrix
  tvec<-t(con)%*%mvec
  tvec<-tvec[,1]
  tempcen<-apply(bcon,1,mean)
  vecz<-rep(0,ncol(con))
  bcon<-t(bcon)
  temp=bcon
  for(ib in 1:nrow(temp))temp[ib,]=temp[ib,]-tempcen+tvec
  smat<-var(temp)
  #smat<-var(bcon-tempcen+tvec)
  chkrank<-qr(smat)$rank
  bcon<-rbind(bcon,vecz)
  if(chkrank==ncol(smat))dv<-mahalanobis(bcon,tvec,smat)
  if(chkrank<ncol(smat)){
    smat<-ginv(smat)
    dv<-mahalanobis(bcon,tvec,smat,inverted=T)
  }
  bplus<-nboot+1
  sig.level<-1-sum(dv[bplus]>=dv[1:nboot])/nboot
  result <- list(test = tvec, p.value = sig.level, call = cl)
  class(result) <- c("spp")
  result
}