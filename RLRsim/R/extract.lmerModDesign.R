#' @importFrom stats model.response
#' @importFrom lme4 getME VarCorr
extract.lmerModDesign <- function(m) {
  X <- getME(m,"X")
  Z <- as.matrix(getME(m,"Z"))
  v <- VarCorr(m)
  resvar <- sigma(m)^2
  Sigma.l <- lapply(v,function(x) x/resvar) #Cov(b)/ Var(Error)
  k <- getME(m,"n_rtrms") #how many grouping factors
  q <- lapply(Sigma.l,NROW) #how many variance components in each grouping factor
  ## OR lapply(m@cnms,length) -- but we should have an extractor for this
  nlevel<-sapply(m@flist, function(x) length(levels(x))) #how many inner blocks in Sigma_i
  ## works as is -- but we should have an extractor
  Vr <- matrix(0,NCOL(Z),NCOL(Z)) #Cov(RanEf)/Var(Error)
  from <- 1
  for(i in 1:k)
  {
    ii<-nlevel[i]
    inner.block<-as.matrix(Sigma.l[[i]])
    to<-from-1+ii*NCOL(inner.block)
    Vr[from:to,from:to]<- inner.block %x% diag(ii)
    from<-to+1
  }
  return(list(
    Vr=Vr, #Cov(RanEf)/Var(Error)
    X=X,
    Z=Z,
    sigmasq=resvar,
    lambda=unique(diag(Vr)),
    y=model.response(model.frame(m)),
    k=k
  ))
}

