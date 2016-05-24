sppbb <- function(formula, id, data, est = "mom", nboot = 500){
#
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
  # Now stack the data in an N by K matrix
  #
  x<-matrix(NA,nrow=nvec[1],ncol=K)
  #
  for(k in 1:K)x[,k]<-xx[[k]]
  kc<-K
  for(j in 2:J){
    temp<-matrix(NA,nrow=nvec[j],ncol=K)
    for(k in 1:K){
      kc<-kc+1
      temp[,k]<-xx[[kc]]
    }
    x<-rbind(x,temp)
  }
  # Now call function rmdzero to do the analysis
  temp<-rmdzero(x,est=est,nboot=nboot)
  result <- list(test = temp$center, p.value = temp$p.value, call = cl)
  class(result) <- c("spp")
  result
}
