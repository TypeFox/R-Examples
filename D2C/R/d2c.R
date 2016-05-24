#' @import MASS randomForest corpcor  lazy


npred<-function(X,Y,lin=TRUE){
  ## normalized mean squared error of the dependency
  N<-NROW(X)
  n<-NCOL(X)
 
  if (n>1){
    w.const<-which(apply(X,2,sd)<0.01)    
    if (length(w.const)>0){
      X<-X[,-w.const]
    }
    n<-NCOL(X)
  }
  
  if (lin)
    return(max(1e-3,regrlin(X,Y)$MSE.loo/var(Y)))
  X<-scale(X)
  e<-Y-lazy.pred(X,Y,X,conPar=c(5,10),
                 linPar=NULL,class=FALSE,cmbPar=10)
  nmse<-mean(e^2)/var(Y) 
  max(1e-3,nmse)
}

#' compute descriptor
#' @param D :  the observed data matrix of size [N,n], where N is the number of samples and n is the number of nodes
#' @param ca : node index (\eqn{1 \le ca \le n}) of the putative cause
#' @param ef : node index (\eqn{1 \le ef \le n}) of the putative effect
#' @param ns : size of the Markov Blanket
#' @param lin : TRUE OR FALSE. if TRUE it uses a linear model to assess a dependency, otherwise a local learning algorithm 
#' @param acc : TRUE OR FALSE. if TRUE it uses the accuracy of the regression as a descriptor
#' @param struct :   TRUE or FALSE to use the ranking in the markov blanket as a descriptor
#' @param pq :  a vector of quantiles used to compute de descriptor
#' @param bivariate :  TRUE OR FALSE. if TRUE it includes the descriptors of the bivariate dependency
#' @details This function is the core of the D2C algorithm. Given two candidate nodes, (\code{ca}, putative cause and \code{ef}, putative effect) it first infers from the dataset D the Markov Blankets of the variables indexed by \code{ca} and \code{ef} (\code{MBca} and \code{MBef}) by using the \link{mimr} algorithm (Bontempi, Meyer, ICML10). Then it computes a set of (conditional) mutual information terms describing the dependency between the variables ca and ef. These terms are used to create a vector of descriptors. If \code{acc=TRUE}, the vector contains the descriptors related to the asymmetric information theoretic terms described in the paper. If \code{struct=TRUE}, the vector contains descriptors related to the positions of the terms of the MBef in MBca and viceversa. The estimation of the information theoretic terms require the estimation of the dependency between nodes. If \code{lin=TRUE} a linear assumption is made. Otherwise the local learning estimator, implemented by the R package \link{lazy}, is used.
#' @references Gianluca Bontempi, Maxime Flauder (2014) From dependency to causality: a machine learning approach. Under submission
#' @references Bontempi G., Meyer P.E. (2010) Causal filter selection in microarray data. ICML'10
#' @references M. Birattari, G. Bontempi, and H. Bersini (1999) Lazy learning meets the recursive least squares algorithm. Advances in Neural Information Processing Systems 11, pp. 375-381. MIT Press.
#' @references G. Bontempi, M. Birattari, and H. Bersini (1999) Lazy learning for modeling and control design. International Journal of Control, 72(7/8), pp. 643-658.
#' @export 
descriptor<-function(D,ca,ef,ns=min(4,NCOL(D)-2),
                     lin=FALSE,acc=TRUE,struct=TRUE, 
                     pq= c(0.1,0.25,0.5,0.75,0.9),
                     bivariate=FALSE ){
  if (bivariate)
    return(c(D2C.n(D,ca,ef,ns,lin,acc,struct,pq=pq),D2C.2(D[,ca],D[,ef])))
  else
    return(c(NCOL(D),NROW(D),D2C.n(D,ca,ef,ns,lin,acc,struct,pq=pq)))
}





norminf<-function(y,x1,x2,lin=TRUE){
  ## I(x1;y| x2)= (H(y|x1)-H(y | x1,x2))/H(y|x1)
  
  np<-npred(x1,y,lin=lin)
  x1x2<-cbind(x1,x2)
  delta<- max(0,np-npred(x1x2,y,lin=lin))/(np+0.01)
  return(delta)
  
}

D2C.n<-function(D,ca,ef,ns=min(4,NCOL(D)-2),
                lin=FALSE,acc=TRUE,struct=TRUE,
                pq= c(0.1,0.25,0.5,0.75,0.9)){
  ## is i cause oj j
  n<-NCOL(D)
  N<-NROW(D)
  
  #### creation of the Markov Blanket of ca (denoted MBca)
  #### MB is obtained by first ranking the other nodes and then selecting a subset of size ns 
  #### with the mimr algorithm
  ind<-setdiff(1:n,ca)
  
  ind<-ind[rankrho(D[,ind],D[,ca],nmax=min(length(ind),50))]
  
  MBca<-ind[mimr(D[,ind],D[,ca],nmax=ns,init=TRUE)]
  
  #### creation of the Markov Blanket of ef (denoted MBef)
  ind<-setdiff(1:n,ef)
  ind<-ind[rankrho(D[,ind],D[,ef],nmax=min(length(ind),50))]
  
  MBef<-ind[mimr(D[,ind],D[,ef],nmax=ns,init=TRUE)]
  
  
  
  namesx<-NULL 
  x<-NULL
  
  
  
  if (struct){
    
    ## position of effect in the MBca 
    if (is.element(ef,MBca))
      pos.ef<-(which(MBca==ef))/(ns)
    else
      pos.ef<-2
    
    ## position of ca in the MBef
    if (is.element(ca,MBef))
      pos.ca<-(which(MBef==ca))/(ns)
    else
      pos.ca<-2
    
    sx.ef<-NULL
    ## position of variables of MBef in MBca 
    for (i in 1:length(MBef)){
      if (is.element(MBef[i],MBca))
        sx.ef<-c(sx.ef,(which(MBca==MBef[i]))/(ns))
      else
        sx.ef<-c(sx.ef,2)
      
    }
    
    ## position of variables of MBca in MBef
    sx.ca<-NULL
    for (i in 1:length(MBca)){
      if (is.element(MBca[i],MBef))
        sx.ca<-c(sx.ca,(which(MBef==MBca[i]))/(ns))
      else
        sx.ca<-c(sx.ca,2)      
    }
    
    
    x<-c(x,pos.ca,pos.ef,quantile(sx.ca,probs=pq),quantile(sx.ef,probs=pq))
    namesx<-c(namesx,"pos.ca","pos.ef",paste0("sx.ca",1:length(pq)),
              paste0("sx.ef",1:length(pq)))
  }
  
  MBca<-setdiff(MBca,ef)
  MBef<-setdiff(MBef,ca)
  
  if (acc){
    
    ## relevance of ca for ef
    ca.ef<-npred(D[,ca],D[,ef],lin=lin)
    ## relevance of ef for ca
    ef.ca<-npred(D[,ef],D[,ca],lin=lin)
    
    
    ## relevance of ca for ef given MBef
    np<-npred(D[,MBef],D[,ef],lin=lin)

    delta<- (npred(D[,c(MBef,ca)],D[,ef],lin=lin)-np)/np
    
    
    ## relevance of ef for ca given MBca
    np<-npred(D[,MBca],D[,ca],lin=lin)
    delta2<- (npred(D[,c(MBca,ef)],D[,ca],lin=lin)-np)/np
    
    
    I1.i<-NULL
    ## Information of Mbef on ca 
    
    for (j in 1:length(MBef)){
      I1.i<-c(I1.i, (npred(D[,MBef[j]],D[,ca],lin=lin)))
    }
    
    I1.j<-NULL
    ## Information of Mbca on ef 
    
    for (j in 1:length(MBca)){
      I1.j<-c(I1.j, (npred(D[,MBca[j]],D[,ef],lin=lin)))
    }
    
    I2.i<-NULL
    ## Information of Mbef on ca given ef
    for (j in 1:length(MBef)){
      I2.i<-c(I2.i, norminf(D[,ca], D[,MBef[j]],D[,ef],lin=lin))
    }
    
    I2.j<-NULL
    ## Information of Mbca on ef given ca
    for (j in 1:length(MBca)){
      I2.j<-c(I2.j, norminf(D[,ef], D[,MBca[j]],D[,ca],lin=lin))
    }
    
    
    I3.i<-NULL
    ## Information of MBef on MBca given ca
    for (i in 1:length(MBca))
      for (j in 1:length(MBef)){
        I3.i<-c(I3.i,(norminf(D[,MBca[i]],D[,MBef[j]],D[,ca],lin=lin)))
      }
    
    I3.j<-NULL
    ## Information of MBef on MBca given ef
    for (i in 1:length(MBca))
      for (j in 1:length(MBef)){
        I3.j<-c(I3.j,(norminf(D[,MBca[i]],D[,MBef[j]],D[,ef],lin=lin)))
      }
    
    
    
    x<-c(x,delta,delta2,ca.ef,ef.ca,
         quantile(I1.i,probs=pq,na.rm=TRUE),quantile(I1.j,probs=pq,na.rm=TRUE),
         quantile(I2.i,probs=pq,na.rm=TRUE),quantile(I2.j,probs=pq,na.rm=TRUE),
         quantile(I3.i,probs=pq,na.rm=TRUE),quantile(I3.j,probs=pq,na.rm=TRUE))
    
    namesx<-c(namesx,"delta","delta2","ca.ef","ef.ca",
              paste0("I1.i",1:length(pq)), paste0("I1.j",1:length(pq)),
              paste0("I2.i",1:length(pq)), paste0("I2.j",1:length(pq)),
              paste0("I3.i",1:length(pq)), paste0("I3.j",1:length(pq)))
  }
  
  if (length(names(x))!=length(namesx))
    browser()
  names(x)<-namesx
  x
  
  
}








D2C.2<-function(x,y,sdkmin=0.5,sdkmax=0.5,Ls=1){
  
  
  
  copula<-TRUE
  alpha<-FALSE
  
  origx<-x
  origy<-y
  
  x<-scale(x)
  y<-scale(y)
  
  sdk<-mean(c(sdkmin,sdkmax))
  
  ##############################@
  N<-length(x)
  point.per.breaks<-min(150,round(N/3))
  if (alpha){
    ncp<-max(2,sqrt(N/point.per.breaks))
    ind.cpx<-seq(min(x),1.01*max(x),length.out=ncp+1)
    ind.cpy<-seq(min(y),1.01*max(y),length.out=ncp+1)
    
    cp<-array(NA,c(ncp,ncp))
    for (i in 1:(ncp))
      for (j in 1:(ncp)){
        
        cp[i,j]<-length(which(x<ind.cpx[i+1] & x>=ind.cpx[i] & y<ind.cpy[j+1] & y>=ind.cpy[j]))/N
        
      }
    
    
    
    alpha<-array(NA,c(ncp-1,ncp-1))
    alpha2<-array(NA,c(ncp-1,ncp-1))
    for (i in 1:(ncp-1))
      for (j in 1:(ncp-1)){
        alpha[i,j]<-log(cp[i,j]*cp[i+1,j+1]/(cp[i,j+1]*cp[i+1,j]))
        
      }
    alpha[which(is.nan(alpha))]<-0
    alpha[which(is.infinite(alpha))]<-0
    alpha[which(alpha==0)]<-NA
    qalpha<-quantile(c(alpha),qnt,na.rm=TRUE)
  }
  
  
  ###################### inp: x1, out: y1, err: Ey
  
  kp<-150
  ix<-sort(x,decreasing=FALSE,ind=TRUE)$ix
  x1<-x[ix]
  y1<-y[ix]
  
  Hy<-NULL
  Ey<-NULL
  ypred<-NULL
  for (i in 1:N ){
    ind<-setdiff(max(1,(i-kp)):min(N,(i+kp)),i)
    
    if (Ls>1){
      S<-seq(sdkmin,sdkmax,length.out=Ls)
      el<-NULL
      for (ssdk in S){
        wd<-numeric(N)
        wd[ind]<- exp(-((x1[ind]-x1[i])^2)/(2*ssdk^2))
        if (sum(wd)>0.01){
          wd<-wd/sum(wd)
          py<-sum(wd*y1)
          el<-c(el,abs(py-y1[i]))
        } else {
          el<-c(el,Inf)
        }
      }
      
      ssdk<-S[which.min(el)]
    } else {
      ssdk<-sdk
    }
    wd<-numeric(N)
    wd[ind]<- exp(-((x1[ind]-x1[i])^2)/(2*ssdk^2))
    if (sum(wd)>0){
      wd<-wd/sum(wd)
      py<-sum(wd*y1)
      sdy<-sqrt(sum(wd*((y1-py)^2)))
    } else {
      py<-mean(y1)
      sdy<-1
    }
    Hy<-c(Hy,sd(y1)-sdy)
    Ey<-c(Ey,y1[i]-py)
    ypred<-c(ypred,py)
  }
  
  EEy<-NULL
  HHy<-NULL
  for (i in 1:N ){
    ind<-setdiff(max(1,(i-kp)):min(N,(i+kp)),i)
    wd<-numeric(N)
    wd[ind]<- exp(-((x1[ind]-x1[i])^2)/(2*sdk^2))
    if (sum(wd)>0){
      wd<-wd/sum(wd)
      
      py<-sum(wd*Ey)
      sdy<-sqrt(sum(wd*((Ey-py)^2)))
    } else {
      py<-mean(Ey)
      sdy<-1
    }
    HHy<-c(HHy,sd(Ey)-sdy)
    
    EEy<-c(EEy,Ey[i]-py)
    
  }
  
  
  ##################################@ inp y2: out x2, err: Ex 
  iy<-sort(y,decreasing=FALSE,ind=TRUE)$ix
  x2<-x[iy]
  y2<-y[iy]
  
  Ex<-NULL
  Hx<-NULL
  xpred<-NULL
  for (i in 1:N ){
    ind<-setdiff(max(1,(i-kp)):min(N,(i+kp)),i)
    if (Ls>1){
      S<-seq(sdkmin,sdkmax,length.out=Ls)
      el<-NULL
      for (ssdk in S){
        wd<-numeric(N)
        wd[ind]<- exp(-((y2[ind]-y2[i])^2)/(2*ssdk^2))
        
        if (sum(wd)>0.01){
          wd<-wd/sum(wd)
          px<-sum(wd*x2)
          el<-c(el,abs(px-x2[i]))
        } else {
          el<-c(el,Inf)
        }
      }
      
      ssdk<-S[which.min(el)]
    } else {
      ssdk<-sdk
    }
    
    wd<-numeric(N)
    wd[ind]<- exp(-((y2[ind]-y2[i])^2)/(2*ssdk^2))
    if (sum(wd)>0){
      wd<-wd/sum(wd)
      px<-sum(wd*x2)
      sdx<-sqrt(sum(wd*((x2-px)^2)))
    } else {
      px<-x2[i]
      sdx<-1
    }
    
    Hx<-c(Hx,sd(x2)-sdx)    
    Ex<-c(Ex,x2[i]-px)
    xpred<-c(xpred,px)
    
  }
  
  EEx<-NULL
  HHx<-NULL
  for (i in 1:N ){
    ind<-setdiff(max(1,(i-kp)):min(N,(i+kp)),i)
    wd<-numeric(N)
    wd[ind]<- exp(-((y2[ind]-y2[i])^2)/(2*sdk^2))
    if (sum(wd)>0){
      wd<-wd/sum(wd)
      px<-sum(wd*Ex)
    } else {
      px<-Ex[i]
      sdx<-1
      
    }
    sdx<-sqrt(sum(wd*((Ex-px)^2)))
    HHx<-c(HHx,sd(Ex)-sdx)
    EEx<-c(EEx,Ex[i]-px)
    
  }
  
  mx<-c(median(Hx),max(Hx))
  dx<-c(max(Hx)-min(Hx))
  
  my<-c(median(Hy),max(Hy))
  dy<-c(max(Hy)-min(Hy))
  
  qnt<-seq(0.05,0.99,length=3)
  
  qx<-quantile(origx,qnt)
  qy<-quantile(origy,qnt)
  cop<- NULL
  ecop<-NULL
  qex<-c(quantile(EEx^2,qnt)-quantile(Ex^2,qnt),
         quantile(Ex,qnt),quantile(Ex^2,qnt))
  qey<-c(quantile(EEy^2,qnt)-quantile(Ey^2,qnt),
         quantile(Ey,qnt),quantile(Ey^2,qnt))
  
  
  
  if (copula){
    
    for (ix in quantile(x,qnt)){
      cx<-NULL
      for (iy in quantile(y,qnt)){
        cx<-c(cx,length(which(x<=ix & y <= iy))/length(x))
        
      }
      cop<-c(cop,(cx))
    }
    
    for (ix in quantile(x,qnt)){
      cx<-NULL
      for (iy in quantile(Ey,qnt)){
        cx<-c(cx,length(which(x<=ix & Ey <= iy))/length(x))
        
      }
      ecop<-c(ecop,(cx))
    }
    
  }
  
  lux<-length(unique(x))
  luy<-length(unique(y))
  
  N<-length(x)
  
  asxy<-c(assoc(x,y),abs(cor(x2,y2))-abs(pcor1(x2,y2,Ex)),
          abs(cor(x1,y1))-abs(pcor1(x1,y1,Ey)))
  asey<-assoc(Ex,y2)
  asex<-assoc(Ex,x2)-asey
  aseyy<-assoc(Ey,y1)
  aseyx<-assoc(Ey,x1)-aseyy
  
  if (sd(Ex)>0.01)
    autx<-c(acf(Ex,plot=FALSE)$acf[2:3],pacf(Ex,plot=FALSE)$acf[1:2])
  else
    autx<-rep(1,4)
  
  if (sd(Ey)>0.01)
    auty<-c(acf(Ey,plot=FALSE)$acf[2:3],pacf(Ey,plot=FALSE)$acf[1:2])
  else
    auty<-rep(1,4)
  
  vxy<-varpred(x,y)
  vyx<-varpred(y,x)
  dv<- c(vxy,vxy-vyx)
  
  nX<-c(paste0("qx",1:length(qx)),paste0("qy",1:length(qy)),
        paste0("dv",1:length(dv)),paste0("cop",1:length(cop)),paste0("N",1:length(N)),
        paste0("mx",1:length(mx)),paste0("dx",1:length(dx)),
        paste0("my",1:length(my)),paste0("dy",1:length(dy)),
        paste0("qex",1:length(qex)),paste0("qey",1:length(qey)),
        paste0("lux",1:length(lux)),paste0("luy",1:length(luy)),
        paste0("asxy",1:length(asxy)),paste0("asex",1:length(asex)),
        paste0("asey",1:length(asey)),paste0("aseyx",1:length(aseyx)),
        paste0("aseyy",1:length(aseyy)),
        paste0("autx",1:length(autx)),paste0("auty",1:length(auty)))
  
  rX<-c(qx,qy,dv,cop,N, mx-my,dx-dy,my,dy,qex-qey,qey,lux,
        luy,asxy,asex,asey,aseyx, aseyy,autx,auty)
  
  names(rX)<-nX
  rX
  
}



varpred<-function(x,y,R=50){
  x<-scale(x)
  y<-scale(y)
  xh<-seq(-1.25,1.25,by=0.1)
  N<-length(x)
  P<-NULL
  beta<-NULL
  for (r in 1:R){
    #set.seed(r)
    ##    wr<-runif(N,0,5)
    ##    wr<-wr/sum(wr)
    Ir<-sample(1:N,round(4*N/5)) ##,replace=TRUE,prob=wr)
    xr<-x[Ir]
    yr<-y[Ir]
    
    px<-NULL
    for ( h in 1:length(xh)){
      sx<-sort(abs(xr-xh[h]),decreasing=FALSE,index=TRUE)$ix[1:min(10,length(xr))]
      px<-c(px,mean(yr[sx]))
      
      
    }
    
    P<-cbind(P,px)
    
  }
  
  
  mean(apply(P,1,sd))
  
  
}
