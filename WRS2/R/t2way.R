t2way <- function(formula, data, tr = 0.2){

  if (missing(data)) {
    mf <- model.frame(formula)
  } else {
    mf <- model.frame(formula, data)
  }
  cl <- match.call()
  
  ## sanity check for incomplete design
  if (any(table(mf[,2], mf[,3]) == 0)) stop("Estimation not possible due to incomplete design.")
  
  J <- nlevels(mf[,2])  # FIXME: convert into factor
  K <- nlevels(mf[,3])
  p <- J*K
  grp <- c(1:p)
  lev.col <- 2:3
  var.col <- 1
  if(tr==.5){
   stop("For medians, use med2way if there are no ties.")
  }

  x <- as.matrix(mf)
  temp=selby2(x,lev.col,var.col)
  #selby(x,lev.col[1],var.col)$grpn
  lev1=length(unique(temp$grpn[,1]))
  lev2=length(unique(temp$grpn[,2]))
  gv=apply(temp$grpn,2,rank)
  gvad=10*gv[,1]+gv[,2]
  grp=rank(gvad)
  J=lev1
  K=lev2
  x=temp$x
  x <- lapply(x, as.numeric)

  tmeans<-0
  h<-0
  v<-0
  for (i in 1:p){
    x[[grp[i]]]=elimna(x[[grp[i]]])
    tmeans[i]<-mean(x[[grp[i]]],tr)
    h[i]<-length(x[[grp[i]]])-2*floor(tr*length(x[[grp[i]]]))
    v[i]<-(length(x[[grp[i]]])-1)*winvar(x[[grp[i]]],tr)/(h[i]*(h[i]-1))
  }
  v<-diag(v,p,p)   # Put squared standard errors in a diag matrix.
  ij<-matrix(c(rep(1,J)),1,J)
  ik<-matrix(c(rep(1,K)),1,K)
  jm1<-J-1
  cj<-diag(1,jm1,J)
  for (i in 1:jm1)cj[i,i+1]<-0-1
  km1<-K-1
  ck<-diag(1,km1,K)
  for (i in 1:km1)ck[i,i+1]<-0-1
  cmat<-kron(cj,ik)  # Contrast matrix for factor A
  alval<-c(1:999)/1000
  for(i in 1:999){
    irem<-i
    Qa<-johan(cmat,tmeans,v,h,alval[i])
    if(Qa$teststat>Qa$crit)break
  }
  A.p.value=irem/1000
  cmat<-kron(ij,ck)  # Contrast matrix for factor B
  for(i in 1:999){
    irem<-i
    Qb<-johan(cmat,tmeans,v,h,alval[i])
    if(Qb$teststat>Qb$crit)break
  }
  B.p.value=irem/1000
  cmat<-kron(cj,ck)  # Contrast matrix for factor A by B
  for(i in 1:999){
   irem<-i
   Qab<-johan(cmat,tmeans,v,h,alval[i])
   if(Qab$teststat>Qab$crit)break
  }
  AB.p.value=irem/1000
  tmeans=matrix(tmeans,J,K,byrow=T)
  result <- list(Qa=Qa$teststat, A.p.value=A.p.value, Qb=Qb$teststat, B.p.value=B.p.value, Qab=Qab$teststat, AB.p.value=AB.p.value, call = cl, varnames = colnames(mf))
  class(result) <- c("t2way")
  result
}
