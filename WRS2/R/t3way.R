t3way <- function(formula, data, tr = 0.2){
  
  if (missing(data)) {
    mf <- model.frame(formula)
  } else {
    mf <- model.frame(formula, data)
  }
  mcl <- match.call()

  J <- nlevels(mf[,2])
  K <- nlevels(mf[,3])
  L <- nlevels(mf[,4])
  p <- J*K*L
  grp <- c(1:p)
  lev.col <- 2:4
  var.col <- 1
  alpha <- 0.05
  
    data=as.matrix(mf)
    temp=selby2(data,lev.col,var.col)
    lev1=length(unique(temp$grpn[,1]))
    lev2=length(unique(temp$grpn[,2]))
    lev3=length(unique(temp$grpn[,3]))
    gv=apply(temp$grpn,2,rank)
    gvad=100*gv[,1]+10*gv[,2]+gv[,3]
    grp=rank(gvad)
    J=lev1
    K=lev2
    L=lev3
    data=temp$x
  
  if(is.matrix(data))data=listm(data)
  data <- lapply(data, as.numeric)
  tmeans<-0
  h<-0
  v<-0
  for (i in 1:p){
    tmeans[i]<-mean(data[[grp[i]]],tr)
    h[i]<-length(data[[grp[i]]])-2*floor(tr*length(data[[grp[i]]]))
    v[i]<-(length(data[[grp[i]]])-1)*winvar(data[[grp[i]]],tr)/(h[i]*(h[i]-1))
  }
  v<-diag(v,p,p)   # Put squared standard errors in a diag matrix.
  ij<-matrix(c(rep(1,J)),1,J)
  ik<-matrix(c(rep(1,K)),1,K)
  il<-matrix(c(rep(1,L)),1,L)
  jm1<-J-1
  cj<-diag(1,jm1,J)
  for (i in 1:jm1)cj[i,i+1]<-0-1
  km1<-K-1
  ck<-diag(1,km1,K)
  for (i in 1:km1)ck[i,i+1]<-0-1
  lm1<-L-1
  cl<-diag(1,lm1,L)
  for (i in 1:lm1)cl[i,i+1]<-0-1
  alval<-c(1:999)/1000
  #  Do test for factor A
  cmat<-kron(cj,kron(ik,il))  # Contrast matrix for factor A
  Qa<-johan(cmat,tmeans,v,h,alpha)
  A.p.value=t3pval(cmat,tmeans,v,h)
  # Do test for factor B
  cmat<-kron(ij,kron(ck,il))  # Contrast matrix for factor B
  Qb<-johan(cmat,tmeans,v,h,alpha)
  B.p.value=t3pval(cmat,tmeans,v,h)
  # Do test for factor C
  cmat<-kron(ij,kron(ik,cl))  # Contrast matrix for factor C
  #Qc<-johan(cmat,tmeans,v,h,alpha)
  for(i in 1:999){
    irem<-i
    Qc<-johan(cmat,tmeans,v,h,alval[i])
    if(Qc$teststat>Qc$crit)break
  }
  C.p.value=irem/1000
  # Do test for factor A by B interaction
  cmat<-kron(cj,kron(ck,il))  # Contrast matrix for factor A by B
  for(i in 1:999){
    irem<-i
    Qab<-johan(cmat,tmeans,v,h,alval[i])
    if(Qab$teststat>Qab$crit)break
  }
  AB.p.value=irem/1000
  # Do test for factor A by C interaction
  cmat<-kron(cj,kron(ik,cl))  # Contrast matrix for factor A by C
  for(i in 1:999){
    irem<-i
    Qac<-johan(cmat,tmeans,v,h,alval[i])
    if(Qac$teststat>Qac$crit)break
  }
  AC.p.value=irem/1000
  #Qac<-johan(cmat,tmeans,v,h,alpha)
  # Do test for factor B by C interaction
  cmat<-kron(ij,kron(ck,cl))  # Contrast matrix for factor B by C
  #Qbc<-johan(cmat,tmeans,v,h,alpha)
  for(i in 1:999){
    irem<-i
    Qbc<-johan(cmat,tmeans,v,h,alval[i])
    if(Qbc$teststat>Qbc$crit)break
  }
  BC.p.value=irem/1000
  # Do test for factor A by B by C interaction
  cmat<-kron(cj,kron(ck,cl))  # Contrast matrix for factor A by B by C
  #Qabc<-johan(cmat,tmeans,v,h,alpha)
  for(i in 1:999){
    irem<-i
    Qabc<-johan(cmat,tmeans,v,h,alval[i])
    if(Qabc$teststat>Qabc$crit)break
  }
  ABC.p.value=irem/1000
  result <- list(Qa=Qa$teststat,A.p.value=A.p.value,
       Qb=Qb$teststat,B.p.value=B.p.value,
       Qc=Qc$teststat,C.p.value=C.p.value,
       Qab=Qab$teststat,AB.p.value=AB.p.value,
       Qac=Qac$teststat,AC.p.value=AC.p.value,
       Qbc=Qbc$teststat,BC.p.value=BC.p.value,
       Qabc=Qabc$teststat,ABC.p.value=ABC.p.value, 
       call = mcl, varnames = colnames(mf))
  class(result) <- c("t3way")
  result
}