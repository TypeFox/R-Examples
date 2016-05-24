gwas = function(y,gen,fam=NULL,chr=NULL,window=NULL,fixed=FALSE){

##################
## INTRODUCTION ##
##################

# REMOVAL OF MISSING Y's  
anyNA = function(x) any(is.na(x))
if(any(is.na(y))){
  w=which(is.na(y))
  y=y[-w];gen=gen[-w,];fam=fam[-w]
}

if(is.null(window)){method="RH"}else{method="EB"}
SNPs=colnames(gen)

# SETTING THE FIXED EFFECT, CHROMOSOME AND FAMILY WHEN IT IS NULL
covariate=matrix(1,length(y),1)
if(is.null(fam)){fam=rep(1,length(y))} else{if(any(is.na(fam))) stop("Family vector must have no NA's")} # calc
if(is.null(chr)){chr=ncol(gen)} # calc

# ORDERING DATA BY FAMILY

cat("Ordering Data",'\n')
organ=function(fam,y,covariate,Z){
  a=cbind(fam,y,covariate,Z);a=a[order(a[,1]),]
  b=array(summary(factor(fam)))
  c=c();for(i in 1:length(b)){d=rep(i,b[i]);c=c(c,d)}
  y=a[,2];fam=c;covariate=a[,3];Z=a[,-c(1:3)];return(list(fam,y,covariate,Z))}
f=organ(fam,y,covariate,gen);fam=f[[1]];y=f[[2]];covariate=f[[3]];gen=f[[4]];rm(f);
covariate=matrix(covariate,ncol=1)

# IMPUTATION OF MISSING Y's - OBSOLETE
  Ymc=function(y,fam){
    MC=function(y){
      x=mean(y,na.rm=T)
      s2=var(y,na.rm=T)
      W=which(is.na(y))
      IMP=rnorm(length(W),x,s2)
      y[W]=IMP;return(y)}
    imp=tapply(y,fam,MC)
    y=unlist(imp);return(as.vector(y))}
  if(any(is.na(y))) y=Ymc(y,fam) # calc

# MARKER IMPUT FUNCTION
gen[gen==5]=NA
IMPUT=function(X){
  X[is.na(X)]=5
  doEXP=function(X){
    Mode=function(x){x=na.omit(x);ux=unique(x);
    ux[which.max(tabulate(match(x,ux)))]}
    exp=Mode(X[,1])
    return(exp)}
  EXP=doEXP(X)
  N=length(X[1,])
  IMP=t(apply(X, MARGIN=1, FUN=inputRow, n = N, exp = EXP))
  return(IMP)}
if(any(is.na(gen))) gen=IMPUT(gen) # calc

# SPARSE DESING MATRIX FUNCTION

cat("Generating Design Matrix",'\n')
SDM = function(gen,fam){
  m = ncol(gen)
  n = nrow(gen)
  u = unique(fam)
  f = length(u)+1
  gg=matrix(0,n,m*f)
  # filling standard parent
  gg[,(1:m-1)*f+1]=gen
  # filling founder parent
  gen=-(gen-2)
  pb=txtProgressBar(style=3)
  for(i in 1:(f-1)){
    w = which(fam==i)
    gg[w,(1:m-1)*f+1+i]=gen[w,]
    setTxtProgressBar(pb,i/(f-1))
  };close(pb)
  return(t(gg))}
GEN=SDM(gen,fam) # calc

# GENOMIC RELATIONSHIP MATRIX
Gmat=function(gen){
  ZZ=crossprod(gen) #individuals are columns
  lambda=sum(diag(ZZ))/nrow(ZZ)
  G=ZZ/lambda
  return(list(G,lambda))}

INPUT=function(gen,fam,chr){
  if(sum(chr)!=ncol(gen)) stop("Number of markers and chromosomes don't match")
  if(length(fam)!=nrow(gen)) stop("Family and genotype don't match")
  MAP=function(gen,fam){ # MGD = Map of Genetic Distance
    ORDER=function(gen,fam){
      Matrix=cbind(fam,gen)
      Matrix=SORT(Matrix)
      mat=Matrix[,-1]
      return(mat)}
    SORT=function(foo){(foo[order(foo[,1]),])}
    # Orgenazing matrices
    Mat=data.frame(ORDER(gen,fam))
    Fam=SORT(matrix(fam,ncol=1))
    nF=dim(array((summary(factor(Fam)))))
    # Functions
    GD=function(snpA,snpB){ # genetic distance
      kosambi = function(r) min(.25*log((1+2*r)/(1-2*r)),.5)
      a0=which(snpA==0)
      a2=which(snpA==2)
      b0=which(snpB==0)
      b2=which(snpB==2)
      NR=length(intersect(a0,b0))+length(intersect(a2,b2))
      RE=length(intersect(a0,b2))+length(intersect(a2,b0))
      if(RE<NR){r=RE/(NR+RE)}else{r=NR/(NR+RE)}
      r = r/(2*(1-r)) # relation r=f(R) in RILs
      k = kosambi(r)
      return(k)}
    AR=function(SNP){ # able to recombine
      ar=function(snp){
        uni=unique(snp)
        int=intersect(uni,c(0,2))
        two=length(int)
        if(two<2) result=NA else result=TRUE
        return(result)}
      ar2=apply(SNP,MARGIN=2,FUN=ar)
      return(ar2)}
    # FUNCTION TO MAP DISTANCE BEFORE SPLITTING
    DOIT=function(Gen){
      nSNP=ncol(Gen)
      GDst=rep(0,nSNP)
      for(i in 2:nSNP){GDst[i]=GD(Gen[,(i-1)],Gen[,i])}
      Good=AR(Gen)
      Good2=c(1,Good[1:(length(Good)-1)])
      GDst=GDst*Good*Good2
      return(GDst)}
    # FUNCTION DOIT BY FAMILY
    DBF=function(gen,fam){ # DOIT by Family
      nF=dim(array((summary(factor(fam)))))
      LIST=split(gen,fam)
      LIST2=lapply(LIST,DOIT)
      A=matrix(0,nF,ncol(gen))
      for(i in 1:nF){A[i,]=unlist(LIST2[i])}
      Final=colMeans(A,na.rm=T)
      return(Final)}
    ccM=function(DBFoutput){
      vect=as.vector(DBFoutput)
      nSNP=length(vect)
      CGD=rep(0,nSNP) # Cumulative gen. dist.
      for(i in 2:nSNP){CGD[i]=vect[i]+CGD[i-1]}
      return(CGD)}
    A=DBF(Mat,Fam)
      A[which(is.na(A))]=0
      for(i in 1:length(A)){if(A[i]>0.5){A[i]=A[i]-0.5}}
    B=ccM(A)
    return(B)}
  # CHROMOSOME VECTOR
  CHROM=function(vector){ 
    len=length(vector);total=sum(vector);
    Vec1=c();
    for(i in 1:len){a=rep(i,vector[i]);Vec1=c(Vec1,a)};
    Vec2=c()
    for(i in 1:len){b=(1:vector[i]);Vec2=c(Vec2,b)};
    Final=cbind(Vec1,Vec2)
    colnames(Final)=c("chr","bin")
    return(Final)}
  a=CHROM(chr)
  ccM=round(MAP(gen,fam),3)
  begin=tapply(ccM,a[,1],function(x){X=x-x[1];return(X)},simplify = T)
  XX=c(); for(i in 1:length(chr)){XX=c(XX,unlist(begin[i]))}; begin=XX; rm(XX)
  end=tapply(begin,a[,1],function(x){X=x[-1];
  X=c(X,(X[length(X)]+0.5));return(X)},simplify = T)
  XX=c(); for(i in 1:length(chr)){XX=c(XX,unlist(end[i]))}; end=XX; rm(XX)
  size=end-begin
  final=cbind(a,begin,end,size,ccM)
  return(final)}
MAP=INPUT(gen,fam,chr) # calc

# Shizhong's MIXED function
mixed<-function(x,y,kk){
  loglike<-function(theta){
    lambda<-exp(theta)
    logdt<-sum(log(lambda*delta+1))
    h<-1/(lambda*delta+1)
    yy<-sum(yu*h*yu)
    yx<-matrix(0,q,1)
    xx<-matrix(0,q,q)
    for(i in 1:q){
      yx[i]<-sum(yu*h*xu[,i])
      for(j in 1:q){
        xx[i,j]<-sum(xu[,i]*h*xu[,j])}}
    loglike<- -0.5*logdt-0.5*(n-q)*log(yy-t(yx)%*%solve(xx)%*%yx)-0.5*log(det(xx))
    return(-loglike)}
  fixed<-function(lambda){
    h<-1/(lambda*delta+1)
    yy<-sum(yu*h*yu)
    yx=timesVec(yu,h,xu,q)
    xx=timesMatrix(xu,h,xu,q,q)
    beta<-solve(xx,yx)
    sigma2<-(yy-t(yx)%*%solve(xx)%*%yx)/(n-q)
    var<-diag(solve(xx)*sigma2)
    stderr<-sqrt(var)
    return(c(beta,stderr,sigma2))}
  qq<-eigen(as.matrix(kk),symmetric=T)
  delta<-qq[[1]]
  uu<-qq[[2]]
  q<-ncol(x)
  n<-ncol(kk)
  vp<-var(y)
  yu<-t(uu)%*%y
  xu<-t(uu)%*%x
  theta<-0
  parm<-optim(par=theta,fn=loglike,hessian = TRUE,method="L-BFGS-B",lower=-5,upper=5)
  lambda<-exp(parm$par)
  conv<-parm$convergence
  fn1<-parm$value
  fn0<-loglike(-Inf)
  lrt<-2*(fn0-fn1)
  hess<-parm$hessian
  parmfix<-fixed(lambda)
  beta<-parmfix[1:q]
  stderr<-parmfix[(q+1):(2*q)]
  sigma2<-parmfix[2*q+1]
  lod<-lrt/4.61
  p_value<-1-pchisq(lrt,1)
  sigma2g<-lambda*sigma2
  goodness<-(vp-sigma2)/vp
  par<-data.frame(conv,fn1,fn0,lrt,goodness,beta,stderr,sigma2,lambda,sigma2g,lod,p_value)
  return(par)}

# Shizhong's BLUP function
blup<-function(gen,map,fam,x,y,kk,beta,lambda,cc){
  qq<-eigen(as.matrix(kk),symmetric=T)
  delta<-qq[[1]]
  uu<-qq[[2]]
  yu<-t(uu)%*%y
  xu<-t(uu)%*%x
  h<-1/(delta*lambda+1)
  r<- max(fam)+1
  m<-nrow(map)
  rr<-yu-xu%*%beta
  gamma<-matrix(0,m,r)
  for(k in 1:m){
    sub<-seq(((k-1)*r+1),((k-1)*r+r))
    z<-as.matrix(gen[sub,])
    zu<-t(uu)%*%t(z)
    zr<-matrix(0,r,1)
    for(i in 1:r){zr[i,1]=sum(rr*h*zu[,i])}
    gamma[k,]<-t(lambda*zr)/cc}
  return(gamma)}

#############
## METHODS ##
#############

RANDOMsma = function (GEN=GEN,MAP=MAP,fam=fam,chr=chr,y=y,COV=covariate,SNPnames=SNPs){
  if(is.vector(y)!=TRUE) stop("Phenotypes: 'y' must be a vector")
  if(is.vector(chr)!=TRUE) stop("Chromosome: 'chr' must be a vector")
  if(is.vector(fam)!=TRUE) stop("Family: 'fam' must be a vector")
  if(is.matrix(gen)!=TRUE) stop("Genotypes: 'gen' must be a matrix")
  if(sum(chr)!=ncol(gen)) stop("Number of markers and chromosomes do not match")
  if(length(fam)!=nrow(gen)) stop("Family and genotype don't match")
  if(length(y)!=length(fam)) stop("Family and phenotypes must have the same length")
  if(length(y)!=nrow(gen)) stop("Dimensions of genotypes and phenotypes do not match")
  map=MAP
  gen=GEN
  cat("Calculating G matrix",'\n')
  G=Gmat(gen)
  kk=(t(G)[[1]]);cc=(t(G)[[2]])
  m<-nrow(map)
  r<-max(fam)+1
  s<-1
  ccM<-map[,6]
  n<-nrow(kk)
  x<-COV
  cat("Solving polygenic model",'\n')
  parm<-mixed(x,y,kk)
  lambda<-parm$lambda
  beta<-parm$beta
  cat("Starting Eigendecomposition",'\n')
  qq<-eigen(as.matrix(kk),symmetric=T)
  delta<-qq[[1]]
  uu<-qq[[2]]
  h<-1/(delta*lambda+1)
  x<-matrix(1,n,1)
  xu<-t(uu)%*%x
  yu<-t(uu)%*%y
  xx<-matrix(0,s,s)
  for(i in 1:s){for(j in 1:s){xx[i,j]<-sum(xu[,i]*h*xu[,j])}}
  loglike<-function(theta){
    xi<-exp(theta)
    tmp0<-zz*xi+diag(r)
    tmp<-xi*solve(tmp0)
    yHy<-yy-t(zy)%*%tmp%*%zy
    yHx<-yx-zx%*%tmp%*%zy
    xHx<-xx-zx%*%tmp%*%t(zx)
    logdt2<-log(det(tmp0))
    loglike<- -0.5*logdt2-0.5*(n-s)*log(yHy-t(yHx)%*%solve(xHx)%*%yHx)-0.5*log(det(xHx))
    return(-loglike)}
  fixed<-function(xi){
    tmp0<-zz*xi+diag(r)
    tmp<-xi*solve(tmp0)
    yHy<-yy-t(zy)%*%tmp%*%zy
    yHx<-yx-zx%*%tmp%*%zy
    xHx<-xx-zx%*%tmp%*%t(zx)
    zHy<-zy-zz%*%tmp%*%zy
    zHx<-zx-zx%*%tmp%*%zz
    zHz<-zz-zz%*%tmp%*%zz
    beta<-solve(xHx,yHx)
    tmp2<-solve(xHx)
    sigma2<-(yHy-t(yHx)%*%tmp2%*%yHx)/(n-s)
    gamma<-xi*zHy-xi*t(zHx)%*%tmp2%*%yHx
    var<-abs((xi*diag(r)-xi*zHz*xi)*as.numeric(sigma2))
    stderr<-sqrt(diag(var))
    result<-list(gamma,stderr,beta,sigma2)
    return(result)}
  parr<-numeric()
  blupp<-numeric()
  cat("Starting Marker Analysis",'\n')
  pb=txtProgressBar(style=3)
  for(k in 1:m){
    sub<-seq(((k-1)*r+1),((k-1)*r+r))
    genk<-gen[sub,]
    zu<-t(uu)%*%t(genk)
    yy<-sum(yu*h*yu)
    zu=as.matrix(zu)
    zx=timesMatrix(xu,h,zu,s,r)
    yx=timesVec(yu,h,xu,s)
    zy=timesVec(yu,h,zu,r)
    zz=timesMatrix(zu,h,zu,r,r)
    theta<-c(0)
    parm<-optim(par=theta,fn=loglike,hessian = TRUE,method="L-BFGS-B",lower=-5,upper=5)
    xi<-exp(parm$par)
    conv<-parm$convergence
    fn1<-parm$value
    fn0<-loglike(-Inf)
    hess<-parm$hessian
    parmfix<-fixed(xi)
    gamma<-parmfix[[1]][1:r]
    stderr<-parmfix[[2]][1:r]
    beta<-parmfix[[3]]
    sigma2<-parmfix[[4]]
    lam_k<-xi
    tau_k<-lam_k*sigma2
    lrt<-2*(fn0-fn1)
    if(lrt<0){
      lrt = 0
      lod = 0
      pval = 0
    }else{
      lod = lrt/4.61
      pval = round(-log(dchisq(lrt,0.5),base = 10),2)
      if(pval<0) pval = 0
    }
    if(r==2){
      names(gamma) = c("allele.eff.A","allele.eff.a")  
      names(stderr) = c("sd.allele.eff.A","sd.allele.eff.a")
    }else{
      names(gamma) = c("allele.eff.standard",paste("allele.eff.founder.",1:(r-1),sep=""))  
      names(stderr) = c("sd.allele.eff.standard",paste("sd.allele.eff.founder.",1:(r-1),sep=""))
    }
    gamma = t(data.frame(gamma))
    stderr = t(data.frame(stderr))
    sigma2g = sigma2*(lambda+lam_k)
    h2 = sigma2g / (sigma2g+sigma2)
    par<-data.frame(conv,fn1,fn0,lod,pval,lrt,sigma2g,sigma2,h2,lam_k,"var.snp"=tau_k,"intercept"=beta,gamma,stderr)
    
    blup<-c(gamma,stderr)
    parr<-rbind(parr,par)
    blupp<-rbind(blupp,blup)
    setTxtProgressBar(pb,k/m)
    #cat(paste("LRT is",round((k/m)*100,2),"% completed",'\n'))
  };close(pb)
  cat("Calculations were performed successfully",'\n')
  rownames(parr) = paste('SNP',1:m,sep='')
  return(list("PolyTest"=parr,"Method"="Empirical Bayes","MAP"=MAP,"SNPs"=SNPnames))}

RandomCIM = function (GEN=GEN,MAP=MAP,fam=fam,chr=chr,y=y,COV=covariate,WIN=window,SNPnames=SNPs){
  if(is.vector(y)!=TRUE) stop("Phenotypes: 'y' must be a vector")
  if(is.vector(chr)!=TRUE) stop("Chromosome: 'chr' must be a vector")
  if(is.vector(fam)!=TRUE) stop("Family: 'fam' must be a vector")
  if(is.matrix(gen)!=TRUE) stop("Genotypes: 'gen' must be a matrix")
  if(sum(chr)!=ncol(gen)) stop("Number of markers and chromosomes do not match")
  if(length(fam)!=nrow(gen)) stop("Family and genotype don't match")
  if(length(y)!=length(fam)) stop("Family and phenotypes must have the same length")
  if(length(y)!=nrow(gen)) stop("Dimensions of genotypes and phenotypes do not match")
  gen=GEN
  map=MAP
  cat("Calculating G matrix",'\n')
  G=Gmat(gen)
  kk=(t(G)[[1]]);cc=(t(G)[[2]])
  m<-nrow(map)
  r<-max(fam)+1
  s<-1
  ccM<-map[,6]
  n<-nrow(kk)
  x<-COV
  cat("Solving polygenic model",'\n')
  parm<-mixed(x,y,kk)
  lambda<-parm$lambda
  beta<-parm$beta
  cat("Starting Eigendecomposition",'\n')
  qq<-eigen(as.matrix(kk),symmetric=T)
  delta<-qq[[1]]
  uu<-qq[[2]]
  h<-1/(delta*lambda+1)
  xu<-t(uu)%*%x
  yu<-t(uu)%*%y
  yy<-sum(yu*h*yu) 
  yx=timesVec(yu,h,xu,s)
  xx=timesMatrix(xu,h,xu,s,s)
  loglike<-function(theta){
    xi<-exp(theta)
    q<-length(xi)
    psi<-NULL
    for(i in 1:q){psi<-c(psi,rep(xi[i],r))}
    psi<-diag(psi)
    tmp0<-zz%*%psi+diag(r*q)
    tmp<-psi%*%solve(tmp0)
    yHy<-yy-t(zy)%*%tmp%*%zy
    yHx<-yx-zx%*%tmp%*%zy 
    xHx<-xx-zx%*%tmp%*%t(zx)
    logdt2<-log(det(tmp0))
    loglike<- -0.5*logdt2-0.5*(n-s)*log(yHy-t(yHx)%*%solve(xHx)%*%yHx)-0.5*log(det(xHx))
    return(-loglike)}
  psi<-NULL
  fixed<-function(xi){q<-length(xi);
                      for(i in 1:q){psi<-c(psi,rep(xi[i],r))}
                      psi<-diag(psi)
                      tmp0<-zz%*%psi+diag(r*q)
                      tmp<-psi%*%solve(tmp0)
                      yHy<-yy-t(zy)%*%tmp%*%zy
                      yHx<-yx-zx%*%tmp%*%zy
                      xHx<-xx-zx%*%tmp%*%t(zx)
                      zHy<-zy-zz%*%tmp%*%zy
                      zHx<-zx-zx%*%tmp%*%zz
                      zHz<-zz-zz%*%tmp%*%zz
                      beta<-solve(xHx,yHx)
                      tmp2<-solve(xHx)
                      sigma2<-(yHy-t(yHx)%*%tmp2%*%yHx)/(n-s)
                      gamma<-psi%*%zHy-psi%*%t(zHx)%*%tmp2%*%yHx
                      var<-abs((psi%*%diag(r*q)-psi%*%zHz%*%psi)*as.numeric(sigma2))
                      stderr<-sqrt(diag(var))
                      result<-list(gamma,stderr,beta,sigma2)
                      return(result)}
  parr<-numeric()
  blupp<-numeric()
  cat("Starting Marker Analysis",'\n')  
  pb=txtProgressBar(style=3)
  window<-WIN
  for(k in 1:m){
    sub<-seq(((k-1)*r+1),((k-1)*r+r))
    gen2<-gen[sub,]
    c1<-ccM[k]-0.5*window
    c3<-ccM[k]+0.5*window
    k1<-suppressWarnings(max(which(ccM<=c1)))
    k3<-suppressWarnings(min(which(ccM>=c3)))
    if(k1==-Inf){gen1<-gen[sub,sample(1:n)]
    } else {
      sub1<-seq(((k1-1)*r+1),((k1-1)*r+r))
      gen1<-gen[sub1,]}
    if(k3==Inf){gen3<-gen[sub,sample(1:n)]
    } else {
      sub3<-seq(((k3-1)*r+1),((k3-1)*r+r))
      gen3<-gen[sub3,]}
    gen1=as.matrix(gen1);gen2=as.matrix(gen2);gen3=as.matrix(gen3);
    ggg=t(rbind(gen1,gen2,gen3));
    zu<-t(uu)%*%ggg
    r3<-3*r
    zx=timesMatrix(xu,h,zu,s,r3)
    zy=timesVec(yu,h,zu,r3)
    zz=timesMatrix(zu,h,zu,r3,r3)
    theta<-c(0,0,0)
    parm<-optim(par=theta,fn=loglike,hessian = TRUE,method="L-BFGS-B",lower=-5,upper=5)
    xi<-exp(parm$par)
    conv<-parm$convergence
    fn1<-parm$value
    hess<-parm$hessian
    parmfix<-fixed(xi)
    gamma<-parmfix[[1]][(r+1):(2*r)]
    stderr<-parmfix[[2]][(r+1):(2*r)]
    beta<-parmfix[[3]]
    sigma2<-parmfix[[4]]
    lam_k<-xi[2]
    tau_k<-lam_k*sigma2
    zx<-matrix(zx[,-seq(r+1,2*r)],s,2*r)
    zy<-matrix(zy[-seq(r+1,2*r),],2*r,1)
    zz<-zz[-seq(r+1,2*r),-seq(r+1,2*r)]
    theta<-c(0,0)
    parm<-optim(par=theta,fn=loglike,hessian = TRUE,method="L-BFGS-B",lower=-5,upper=5)
    fn0<-parm$value
    lrt<-2*(fn0-fn1)
    if(lrt<0){
      lrt = 0
      lod = 0
      pval = 0
    }else{
      lod = lrt/4.61
      pval = round(-log(dchisq(lrt,0.5),base = 10),2)
      if(pval<0) pval = 0
    }    
    if(r==2){
      names(gamma) = c("allele.eff.A","allele.eff.a")  
      names(stderr) = c("sd.allele.eff.A","sd.allele.eff.a")
    }else{
      names(gamma) = c("allele.eff.standard",paste("allele.eff.founder.",1:(r-1),sep=""))  
      names(stderr) = c("sd.allele.eff.standard",paste("sd.allele.eff.founder.",1:(r-1),sep=""))
    }
    gamma = t(data.frame(gamma))
    stderr = t(data.frame(stderr))
    sigma2g = sigma2*(lambda+lam_k)
    h2 = sigma2g / (sigma2g+sigma2)
    par<-data.frame(conv,fn1,fn0,lod,pval,lrt,sigma2g,sigma2,h2,lam_k,"var.snp"=tau_k,"intercept"=beta,gamma,stderr)
    blup<-c(gamma,stderr)
    parr<-rbind(parr,par)
    blupp<-rbind(blupp,blup)
    setTxtProgressBar(pb,k/m)
    #cat("LRT",round((k/m)*100,2),"% completed",'\n')
  };close(pb)
  cat("Calculations were performed successfully",'\n')
  rownames(parr) = paste('SNP',1:m,sep='')
  return(list("PolyTest"=parr,"Method"="Empirical Bayes with moving-window strategy","MAP"=MAP,"SNPs"=SNPnames))}

FIXEDsma = function (GEN=GEN,MAP=MAP,fam=fam,chr=chr,y=y,COV=covariate,SNPnames=SNPs){
  if(is.vector(y)!=TRUE) stop("Phenotypes: 'y' must be a vector")
  if(is.vector(chr)!=TRUE) stop("Chromosome: 'chr' must be a vector")
  if(is.vector(fam)!=TRUE) stop("Family: 'fam' must be a vector")
  if(is.matrix(gen)!=TRUE) stop("Genotypes: 'gen' must be a matrix")
  if(sum(chr)!=ncol(gen)) stop("Number of markers and chromosomes do not match")
  if(length(fam)!=nrow(gen)) stop("Family and genotype don't match")
  if(length(y)!=length(fam)) stop("Family and phenotypes must have the same length")
  if(length(y)!=nrow(gen)) stop("Dimensions of genotypes and phenotypes do not match")
  map=MAP
  gen=GEN
  cat("Calculating genomic kinship",'\n')
  G=Gmat(gen)
  kk=(t(G)[[1]]);cc=(t(G)[[2]])
  m<-nrow(map)
  r<-max(fam)+1
  q<-1
  ccM<-map[,6]
  n<-nrow(kk)
  x<-COV
  cat("Solving polygenic model (P3D)",'\n')
  parm<-mixed(x,y,kk)
  lambda<-parm$lambda
  beta<-parm$beta
  cat("Starting Eigendecomposition",'\n')
  qq<-eigen(as.matrix(kk),symmetric=T)
  delta<-qq[[1]]
  uu<-qq[[2]]
  h<-1/(delta*lambda+1)
  xu<-t(uu)%*%x
  yu<-t(uu)%*%y
  yy<-sum(yu*h*yu)
  yx = timesVec(yu,h,xu,q)
  xx = timesMatrix(xu,h,xu,q,q)
  beta<-solve(xx,yx)
  parr<-numeric()
  blupp<-numeric()
  cat("Starting Single Marker Analysis",'\n')
  pb=txtProgressBar(style=3)
  for(k in 1:m){
    sub<-seq(((k-1)*r+1),((k-1)*r+r))
    zu<-t(uu)%*%t(gen[sub,])
    zu = as.matrix(zu)
    zx = timesMatrix(xu,h,zu,q,r)
    zy = timesVec(yu,h,zu,r)
    zz = timesMatrix(zu,h,zu,r,r)
    zzi<-solve(zz+diag(.000001,ncol(zz)))
    b<-solve(zz+diag(.000001,ncol(zz)),zy)
    s2<-(yy-t(zy)%*%zzi%*%zy)/(n-r-q);
    v<-zzi*drop(s2);
    g<-b-mean(b)
    gamma<-g
    wald<-t(g)%*%solve(v)%*%g
    stderr<-sqrt(diag(v))
    sigma2<-s2
    if(r==2){
      rownames(gamma) = c("allele.eff.A","allele.eff.a")  
      names(stderr) = c("sd.allele.eff.A","sd.allele.eff.a")
    }else{
      rownames(gamma) = c("allele.eff.standard",paste("allele.eff.founder.",1:(r-1),sep=""))  
      names(stderr) = c("sd.allele.eff.standard",paste("sd.allele.eff.founder.",1:(r-1),sep=""))
    }
    gamma = t(data.frame(gamma))
    stderr = t(data.frame(stderr))
    par<-data.frame(wald,sigma2,gamma,stderr)
    parr<-rbind(parr,par)
    setTxtProgressBar(pb,k/m)
  };close(pb)
  rownames(parr) = paste('SNP',1:m,sep='')
  cat("Calculations were performed successfully",'\n')
  return(list("PolyTest"=parr,"Method"="P3D","MAP"=MAP,"SNPs"=SNPnames))}


#########
## RUN ##
#########

if(fixed==TRUE){
  fit=FIXEDsma(GEN=GEN,MAP=MAP,fam=fam,chr=chr,y=y,COV=covariate,SNPnames=SNPs)
}else{
  if(method=="RH"){fit=RANDOMsma(GEN=GEN,MAP=MAP,fam=fam,chr=chr,y=y,COV=covariate,SNPnames=SNPs)}
  if(method=="EB"){fit=RandomCIM(GEN=GEN,MAP=MAP,fam=fam,chr=chr,y=y,COV=covariate,WIN=window,SNPnames=SNPs)}
  neg = which(fit$PolyTest$lrt<0);fit$PolyTest$lrt[neg]=0
}
class(fit) <- "NAM"

return(fit)}
