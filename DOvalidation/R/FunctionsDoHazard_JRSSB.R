################################################################################### 
#  Functions to perform the calculations in the paper:
#         Do-validating local linear hazards 
#  by Gamiz, Mammen, Martinez-Miranda and Nielsen (2014)
###################################################################################

# Epanechnikov kernel
K.epa<-function(u){return(0.75*(1-u^2)*(abs(u)<1))}
# Sextic kernel
K.sextic<-function(u) { return(((3003/2048)*(1-(u)^2)^6)*(abs(u)<1))}


## Local linear hazard estimator from discrete data 
hazard.LL<-function(xi,Oi,Ei,x,b,K="epa",Ktype="symmetric",CI=FALSE)
{
  #we have discretized data:
  # xi are the grid points (mid points of the considered intervals)
  # Oi is the vector with the number of occurrences into the interval
  # Ei is the vector with the exposure into the interval per unit of time,
  #     this is Yi*(ti-t(i-1)) with Yi being the total number of individual at risk
  #     at the beginning of the interval
  # x is the estimation point (or vector)
  # b is the bandwidth parameter (scalar)
  # K is the kernel function
  # Ktype is one of "symmetric", "left" or "right" so Ku is considered as
  #   the standard symmetric kernel or the one-sided versions
  
  M<-length(xi)
  ngrid<-length(x)
  if (ngrid>1) {u<-sapply(1:M,function(i) (x-xi[i])/b); u<-t(u)} else {u<-(x-xi)/b}
  uu<-u*b 
  if (K=="epa"){
    Ku<-function(u) K.epa(u)
    if (Ktype=="left") Ku<-function(u) K.epa(u)*2*( (abs(u)<1) & (u<0)  )
    if (Ktype=="right") Ku<-function(u) K.epa(u)*2*( (abs(u)<1) & (u>0)  )
  }
  if (K=="sextic"){
    Ku<-function(u) K.sextic(u)
    if (Ktype=="left") Ku<-function(u) K.sextic(u)*2*( (abs(u)<1) & (u<0)  )
    if (Ktype=="right") Ku<-function(u) K.sextic(u)*2*( (abs(u)<1) & (u>0)  )
  }
  if (ngrid>1)
     {
      a1<-b^(-1)* colSums(Ku(u)*uu*Ei,na.rm=TRUE)
      a2<-b^(-1)*colSums(Ku(u)*(uu^2)*Ei,na.rm=TRUE)
      }
  else 
     {
      a1<-b^(-1)*sum(Ku(u)*uu*Ei,na.rm=TRUE)
      a2<-b^(-1)*sum(Ku(u)*(uu^2)*Ei,na.rm=TRUE)
     }
  if (ngrid>1) {correct<-t(apply(uu,1,function(y){return(a2-y*a1)}))}
     else {correct<-(a2-a1*uu)}
  
  Kequiv.u<-correct*Ku(u)

  if (ngrid>1)  OLL<-b^(-1)*colSums(Kequiv.u*Oi,na.rm=T) else {
  	  OLL<-b^(-1)*sum(Kequiv.u*Oi,na.rm=T)} 
  
  if (ngrid>1)  ELL<-b^(-1)*colSums(Kequiv.u*Ei,na.rm=T) else {
  	  ELL<-b^(-1)*sum(Kequiv.u*Ei,na.rm=T)} 

  hLL<-OLL/ELL
  hLL[ELL==0]<-NA

  # normalized ocurrences and exposures    
  if (ngrid>1)  sum.K<-b^(-1)*colSums(Kequiv.u,na.rm=T) else {
    sum.K<-b^(-1)*sum(Kequiv.u,na.rm=T)}  
  OLL.norm<-OLL/sum.K
  ELL.norm<-ELL/sum.K
  
  # Pointwise CI 95%
  if (CI==FALSE) {CI.inf<-NA;CI.sup<-NA} else {
    Yi.smooth<-ELL.norm
    K2<-function(u){return(Ku(u)^2)}
    R.K<-integrate(K2,lower=-1,upper=1)$value
    #The limits of the confidence interval are:
    z1<-qnorm(0.025); z2<-qnorm(0.975)
    limit.inf.ti<-z1*sqrt((R.K*hLL)/(b*Yi.smooth));
    limit.sup.ti<-z2*sqrt((R.K*hLL)/(b*Yi.smooth));
    CI.inf<-hLL+limit.inf.ti
    CI.sup<-hLL+limit.sup.ti
  }

  hazard.res<-list(x=x,OLL=OLL,ELL=ELL,hLL=hLL,OLL.norm=OLL.norm,ELL.norm=ELL.norm,CI.inf=CI.inf,CI.sup=CI.sup)
  
  return(hazard.res)
}

########################################################################
## The cross-validation bandwidth estimate

b.CV<-function(grid.b,nb=50,K="epa",xi,Oi,Ei) 
{
  M<-length(xi)
  if (missing(grid.b)){
    amp<-xi[M]-xi[1]
    b.min<-amp/(M+1)
    b.max<-amp/2
    grid.b<-seq(b.min,b.max,length=nb)
  } else {nb<-length(grid.b)}
  
  cv.score<-function(b)
  {
    alpha.i<-hazard.LL(xi,Oi,Ei,x=xi,b,K,Ktype="symmetric")$hLL
    cv1<-sum((alpha.i^2)*Ei,na.rm=TRUE)
    O.i<-matrix(Oi,M,M); 
    i0<-which(Oi==0)
	  Oi0<-Oi-1; Oi0[i0]<-0
    diag(O.i)<-Oi0 ; 
    hi.xi<-sapply(1:M, function(i) {hazard.LL(xi,O.i[,i],Ei,x=xi[i],b,K,Ktype="symmetric")$hLL})
    cv2<-sum(hi.xi*Oi,na.rm=TRUE)
    cv<-(cv1-2*cv2)
    cv[cv==0]<-NA
    return(cv)  
  }

  cv.values<-sapply(grid.b,cv.score)
  ind.cv<-which.min(cv.values)
  if ((ind.cv==1)|(ind.cv==nb)) warning("The CV score doesn't have a minumum in the grid of bandwidths")
  bcv<-grid.b[ind.cv]
  cv.res<-list(bcv=bcv,ind.cv=ind.cv,cv.values=cv.values,grid.b=grid.b)
  return(cv.res)
}

## The one-sided cross-validation bandwidth estimate

b.OSCV<-function(grid.b,nb=50,K="epa",Ktype="left",xi,Oi,Ei) 
{
  # The rescaling constant for do-validation
  if (K=="epa") {Cval<-0.5232}
  if (K=="sextic") {Cval<-0.5874 }
  M<-length(xi)
  if (missing(grid.b)){
    amp<-xi[M]-xi[1]
    b.min<-amp/(M+1)
    b.max<-amp
    grid.b<-seq(b.min,b.max,length=nb)
  } 
   
  cv.score<-function(b)
  {
    alpha.i<-hazard.LL(xi,Oi,Ei,x=xi,b,K,Ktype)$hLL
    cv1<-sum((alpha.i^2)*Ei,na.rm=TRUE)
    O.i<-matrix(Oi,M,M);
    i0<-which(Oi==0)
    Oi0<-Oi-1; Oi0[i0]<-0
    diag(O.i)<-Oi0  
    hi.xi<-sapply(1:M, function(i) {hazard.LL(xi,O.i[,i],Ei,x=xi[i],b,K,Ktype)$hLL})
    cv2<-sum(hi.xi*Oi,na.rm=TRUE)
    cv<-(cv1-2*cv2)
    cv[cv==0]<-NA
    return(cv)  
  }
  
  oscv.values<-sapply(grid.b,cv.score) 
  
  ind.oscv<-which.min(oscv.values)
  boscv<-grid.b[ind.oscv]*Cval ## Rescale to the original grid of bandwidths
  if ((ind.oscv==1)|(ind.oscv==nb)) warning("The OSCV score doesn't have a minumum in the grid of bandwidths")
  oscv.res<-list(boscv=boscv,ind.oscv=ind.oscv,oscv.values=oscv.values,Cval=Cval,
                 grid.b=grid.b)
  return(oscv.res)
}


## Local linear hazard with RAMLAU-HANSEN weighting

hazard.LL.RH<-function(xi,Oi,Ei,x,b,K="epa")
{  
  M<-length(xi)
  #calculate the mean points in each interval
  deltai<-diff(c(39,xi))
  Yi<-Ei*deltai  
  ngrid<-length(x)
  if (ngrid>1) {u<-sapply(1:M,function(i) (x-xi[i])/b); u<-t(u)} else {
    u<-(x-xi)/b
  }
  uu<-u*b
  
  if (K=="epa") Ku<-function(u) K.epa(u)
  if (K=="sextic") Ku<-function(u) K.sextic(u)
  
  if (ngrid>1)
  {c0.d<-b^(-1)*colSums(Ku(u),na.rm=T) 
   c1.d<-b^(-1)* colSums(Ku(u)*uu,na.rm=T)
   c2.d<-b^(-1)*colSums(Ku(u)*(uu^2),na.rm=T)
  }
  else 
  {c0.d<-b^(-1)*sum(Ku(u)*deltai,na.rm=T) 
   c1.d<-b^(-1)*sum(Ku(u)*uu*deltai,na.rm=T)
   c2.d<-b^(-1)*sum(Ku(u)*(uu^2)*deltai,na.rm=T)
  }
    
  if (ngrid>1) {correct<-t(apply(uu,1,function(y){return((c2.d-y*c1.d)/(c2.d*c0.d-c1.d^2))}))}
  else {correct<-((c2.d-c1.d*uu)/(c2.d*c0.d-c1.d^2))}
  
  Ku.u<-correct*Ku(u)
  hi<-Oi/Yi
  hi[Yi==0]<-NA
  
  if (ngrid>1)  RH.d<-b^(-1)*colSums(Ku.u*hi,na.rm=T) else {
    RH.d<-b^(-1)*sum(Ku.u*hi,na.rm=T)} 

  hazard.res<-list(x=x,hLL=RH.d)
  
  return(hazard.res)
}

