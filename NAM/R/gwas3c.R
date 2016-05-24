gwas2 = function(y,gen,fam=NULL,chr=NULL,fixed=FALSE,EIG=NULL,cov=NULL){
  
  ##################
  ## INTRODUCTION ##
  ##################
  
  # REMOVAL OF MISSING Y's  
  anyNA = function(x) any(is.na(x))
  if(any(is.na(y))){
    wMIS=which(is.na(y))
    y=y[-wMIS];gen=gen[-wMIS,];fam=fam[-wMIS];cov=cov[-wMIS]
  }else{wMIS=NULL}
  
  method="RH"
  fx=fixed
  SNPs=colnames(gen)
  
  # SETTING THE FIXED EFFECT, CHROMOSOME AND FAMILY WHEN IT IS NULL
  if(is.null(cov)){y=y-mean(y);covariate=matrix(1,length(y),1)} else covariate=matrix(cov,ncol=1)
  if(is.null(fam)){fam=rep(1,length(y))} else{if(any(is.na(fam))) stop("Family vector must have no NA's")} # calc
  if(is.null(chr)){chr=ncol(gen)} # calc
  
  # ORDERING DATA BY FAMILY
  if(is.data.frame(gen)) gen=data.matrix(gen)
  
  if(mean(order(fam)==1:length(fam))!=1){
    cat("Ordering Data",'\n')
    W = order(fam)
    fam = fam[W]
    y = y[W]
    covariate = covariate[W]
    gen = gen[W,]
    if(!is.null(EIG)){EIG$vectors = EIG$vectors[W,]}
  }
  covariate=matrix(covariate,ncol=1)  
  
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
  if(any(is.na(gen))){gen=IMPUT(gen)} # calc
  
  # acceleration
  GGG=function(G,fam){
    f=max(fam)+1
    POP = function(gfa){
      J = rep(0,f);g = gfa[1];fa = gfa[2]
      if(g==2){J[1]=2}else{if(g==1)
      {J[1]=1;J[fa+1]=1}else{J[fa+1]=2}}
      return(J)}
    gfa = cbind(G,fam)
    return(unlist(apply(gfa,1,POP)))}
  
  if(any(fam!=1)){
    Gmat = function(gen,fam){
      # common parent linear kernel
      g1 = tcrossprod(gen)
      g1 = g1/mean(diag(g1))  
      # founder parent linear kernel
      g2 = ((gen-1)*-1)+1
      g2 = tcrossprod(g2)
      g2 = g2/mean(diag(g2))
      # adjusting intra-family relationship
      for(i in unique(fam)){
        nam = which(fam==i)
        g1[nam,nam]=g2[nam,nam]+g1[nam,nam]}
      # Final estimates
      lambda = mean(diag(g1))
      G = g1/lambda
      return(list(G,lambda))}
    
  }else{
    Gmat = function(gen,fam){
      # common parent linear kernel
      g1 = tcrossprod(gen)
      # Final estimates
      lambda = mean(diag(g1))
      G = g1/lambda
      return(list(G,lambda))}
  }
  
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
        d = kosambi(r)
        return(d)}
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
    
    if(!is.null(EIG)){
      qq=EIG
      if(!is.null(wMIS)){
        len = length(qq$values)-length(wMIS)
        qq$values = qq$values[1:len]*(len/length(qq$values))
        qq$vectors = qq$vectors[-wMIS,]
        CE = which(order(abs(cor(qq$vectors,y)))<=len)
        qq$vectors = qq$vectors[,CE]
      }
    }else{
      qq=eigen(as.matrix(kk),symmetric=T)
    }
    
    delta<-qq[[1]]
    uu<-qq[[2]]
    q<-ncol(x)
    n<-ncol(uu)
    vp<-var(y)
    yu<-crossprod(uu,y)
    xu<-crossprod(uu,x)
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
    if(!is.null(EIG)){
      qq=EIG
      if(!is.null(wMIS)){
        len = length(qq$values)-length(wMIS)
        qq$values = qq$values[1:len]*(len/length(qq$values))
        qq$vectors = qq$vectors[-wMIS,]
        CE = which(order(abs(cor(qq$vectors,y)))<=len)
        qq$vectors = qq$vectors[,CE]
      }
    }else{qq<-eigen(as.matrix(kk),symmetric=T)}
    delta<-qq[[1]]
    uu<-qq[[2]]
    yu<-crossprod(uu,y)
    xu<-crossprod(uu,x)
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
  
  RANDOMsma = function (GEN=gen,MAP=MAP,fam=fam,chr=chr,y=y,COV=covariate,SNPnames=SNPs){
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
    
    if(is.null(EIG)){
      cat("Calculating genomic kinship\n")
      G=Gmat(gen,fam)
      kk=(t(G)[[1]])
      cc=(t(G)[[2]])
    }else{
      kk=1
      cc=2
    }
    
    m<-nrow(map)
    r<-max(fam)+1
    s<-1
    n<-length(y)
    x<-COV
    cat("Solving polygenic model",'\n')
    parm<-mixed(x,y,kk)
    lambda<-parm$lambda
    beta<-parm$beta
    if(!is.null(EIG)){
      qq=EIG
      if(!is.null(wMIS)){
        len = length(qq$values)-length(wMIS)
        qq$values = qq$values[1:len]*(len/length(qq$values))
        qq$vectors = qq$vectors[-wMIS,]
        CE = which(order(abs(cor(qq$vectors,y)))<=len)
        qq$vectors = qq$vectors[,CE]
      }
    }else{
      cat("Starting Eigendecomposition",'\n')
      qq<-eigen(as.matrix(kk),symmetric=T)}
    delta<-qq[[1]]
    uu<-qq[[2]]
    h<-1/(delta*lambda+1)
    yu<-crossprod(uu,y)
    xu<-crossprod(uu,x)
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
    if(max(fam)==1) r=1
    parr<-matrix(NA,nrow = m,ncol = (12+r))
    cat("Starting Marker Analysis",'\n')
    pb=txtProgressBar(style=3)
    for(k in 1:m){
      if(r==2){
        genk<-t(gen[,k])
      }else{genk<-GGG(gen[,k],fam)}
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
      lrt<- 2*(fn0-fn1)
      if(lrt<0){
        lrt = 0
        lod = 0
        pval = 0
      }else{
        lod = lrt/4.61
        pval = round(-log(dchisq(lrt,0.5),base = 10),2)
        if(pval<0) pval = 0
      }
      if(r>2) names(gamma) = c("allele.eff.standard",paste("allele.eff.founder.",1:(r-1),sep=""))
      gamma = t(data.frame(gamma))
      sigma2g = sigma2*(lambda+lam_k)
      h2 = sigma2g / (sigma2g+sigma2)
      par<-data.frame(conv,fn1,fn0,lod,pval,lrt,sigma2g,sigma2,h2,lam_k,"var.snp"=tau_k,"intercept"=beta,gamma)
      parr[k,] = as.numeric(par)
      setTxtProgressBar(pb,k/m)
    }
    parr=data.frame(parr)
    names(parr) = names(par)
    close(pb)
    cat("Calculations were performed successfully",'\n')
    if(max(fam)==1){names(parr)[13] = c("allele.eff")}
    rownames(parr) = paste('SNP',1:m,sep='')
    return(list("PolyTest"=parr,"Method"="Empirical Bayes","MAP"=MAP,"SNPs"=SNPnames))}
  
  FIXEDsma = function (GEN=gen,MAP=MAP,fam=fam,chr=chr,y=y,COV=covariate,SNPnames=SNPs){
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
    
    if(is.null(EIG)){
      cat("Calculating genomic kinship \n")
      G=Gmat(gen,fam)
      kk=(t(G)[[1]])
      cc=(t(G)[[2]])
    }else{
      kk=1
      cc=2
    }
    
    m<-nrow(map)
    r<-max(fam)+1
    q<-1
    ccM<-map[,6]
    n<-length(y)
    x<-COV
    cat("Solving polygenic model",'\n')
    parm<-mixed(x,y,kk)
    lambda<-parm$lambda
    beta<-parm$beta
    if(!is.null(EIG)){
      qq=EIG
      if(!is.null(wMIS)){
        len = length(qq$values)-length(wMIS)
        qq$values = qq$values[1:len]*(len/length(qq$values))
        qq$vectors = qq$vectors[-wMIS,]
        CE = which(order(abs(cor(qq$vectors,y)))<=len)
        qq$vectors = qq$vectors[,CE]
      }
    }else{
      cat("Starting Eigendecomposition",'\n')
      qq<-eigen(as.matrix(kk),symmetric=T)}
    delta<-qq[[1]]
    uu<-qq[[2]]
    h<-1/(delta*lambda+1)
    yu<-crossprod(uu,y)
    xu<-crossprod(uu,x)
    yy<-sum(yu*h*yu)
    yx = timesVec(yu,h,xu,q)
    xx = timesMatrix(xu,h,xu,q,q)
    beta<-solve(xx,yx)
    if(max(fam)==1) r=1
    parr<-matrix(NA,nrow=m,ncol=(2+r))
    cat("Starting Single Marker Analysis",'\n')
    pb=txtProgressBar(style=3)
    for(k in 1:m){
      if(r==1){
        genk<-t(gen[,k])
      }else{genk<-GGG(gen[,k],fam)}
      zu<-t(uu)%*%t(genk)
      zu = as.matrix(zu)
      zx = timesMatrix(xu,h,zu,q,r)
      zy = timesVec(yu,h,zu,r)
      zz = timesMatrix(zu,h,zu,r,r)
      zzi<-solve(zz+diag(.000001,ncol(zz)))
      b<-solve(zz+diag(.000001,ncol(zz)),zy)
      s2<-(yy-t(zy)%*%zzi%*%zy)/(n-r-q)
      v<-zzi*drop(s2);
      if(r==1){
        g<-b-mean(b); gamma<-g
        wald<-t(b)%*%solve(v)%*%b
      }else{
        g<-b-mean(b); gamma<-g
        wald<-t(g)%*%solve(v)%*%g}
      sigma2<-s2
      if(r>1){
        rownames(gamma) = paste("allele.eff.",1:r,sep="")
        gamma = t(data.frame(gamma))
      }   
      par = data.frame(wald,sigma2,gamma)
      parr[k,] = as.numeric(par)
      setTxtProgressBar(pb,k/m)
    }
    close(pb)
    parr = data.frame(parr)
    names(parr) = names(par)
    if(r==1) names(parr)[3] = c("allele.eff")
    cat("Calculations were performed successfully",'\n')
    rownames(parr) = paste('SNP',1:m,sep='')
    return(list("PolyTest"=parr,"Method"="P3D","MAP"=MAP,"SNPs"=SNPnames))}
  
  #########
  ## RUN ##
  #########
  
  if(fx==1){
    fit=FIXEDsma(GEN=gen,MAP=MAP,fam=fam,chr=chr,y=y,COV=covariate,SNPnames=SNPs)
  }else{
    fit=RANDOMsma(GEN=gen,MAP=MAP,fam=fam,chr=chr,y=y,COV=covariate,SNPnames=SNPs)
    moda=function(x){
      it=5;ny=length(x);k=ceiling(ny/2)-1; while(it>1){
        y=sort(x); inf=y[1:(ny-k)]; sup=y[(k+1):ny]
        diffs=sup-inf; i=min(which(diffs==min(diffs)))
        M=median(y[i:(i+k)]); it=it-1}; return(M)}
    neg = which(fit$PolyTest$lrt<0);fit$PolyTest$lrt[neg]=0
    mo = moda(fit$PolyTest$lrt)
    mos = which(fit$PolyTest$lrt==mo);fit$PolyTest$lrt[mo]=0
  }
  class(fit) <- "NAM"
  return(fit)}

eigX = function(gen,fam){
  if(any(is.na(gen))){
    mis = apply(gen,2,function(x) mean(is.na(x)))
    gen = gen[,mis<.8]
    gen = markov(gen,ncol(gen))
  }
  gm=function(gen,fam){
    g1 = tcrossprod(gen)
    g1 = g1/mean(diag(g1))
    g2 = ((gen-1)*-1)+1
    g2 = tcrossprod(g2)
    g2 = g2/mean(diag(g2))
    for (i in unique(fam)){nam = which(fam == i)
    g1[nam, nam] = g2[nam, nam] + g1[nam, nam]}
    lambda = mean(diag(g1))
    G = g1/mean(diag(g1))
    return(G)}
  return(eigen(gm(gen,fam),symmetric=TRUE))
}