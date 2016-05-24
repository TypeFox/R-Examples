#' @include gridUtilities.r
#' @include tailIntEmp.r
NULL

tailSmithInt <- function(loc, siginv){
    locdiff<-loc[3:4]-loc[1:2]
    a2 <- (t(locdiff) %*% siginv %*% locdiff)
    a <- sqrt(a2)
    return(as.vector(pnorm(a/2) + (exp(a2)*pnorm((-3*a)/2))/3))
}

tominimizeSmith<-function(sigvector,totlist,pairs,w){
  sigm <- rbind(c(sigvector[1], sigvector[3]), c(sigvector[3], sigvector[2]))
  siginv <- ((1/(sigm[1,1]*sigm[2,2] - sigm[1,2]*sigm[2,1]))* matrix(c(sigm[2,2], -sigm[2,1], -sigm[1,2], sigm[1,1]), nrow=2))
  res<-apply(pairs,1,tailSmithInt,siginv=siginv) - totlist
  return(t(res) %*% w %*% res)
}

ellSm <- function(t, loc, siginv){ 
  if((loc[1]==0) && (loc[2]==0)){
    return(pmax(t[1],t[2]))
  } else{
    a<-sqrt((t(loc)%*%siginv%*%loc))
    res<-.C("ellsmith2",as.double(t),as.double(a),result=double(1),PACKAGE="spatialTailDep")$result
    return(res)
  }
}

ellSm4d<-function(t,siglist,a2,d){ 
  result<-vector(length=d)
  for(i in 1:d){
    const<-(a2[[i]] + log(t[i]/t[-i]))
    result[i]<-t[i]*pmvnorm(lower=rep(-Inf,d-1),upper=const,sigma=siglist[[i]],algorithm=TVPACK)[1]
  }
  return(sum(result))
}

tailSmith <- function(t, loc, siginv,d){ 
    if(d==2){
      return(ellSm(t,loc[,2]-loc[,1],siginv))
    } else{
      siglist<-a2<-vector('list',length=d)
      for(i in 1:d){
        siglist[[i]]<-t(loc[,i]%*%t(rep(1,d-1)) - loc[,-i])%*%siginv%*%(loc[,i]%*%t(rep(1,d-1)) - loc[,-i])
        a2[[i]]<-diag(0.5*t(loc[,i]-loc[,-i])%*%siginv%*%(loc[,i]-loc[,-i]))
      }
      return(ellSm4d(t,siglist,a2,d))
    }
}

ellSm4dv2<-function(t,siglist,a2,d,ind){ 
  return(2*t[ind]*ellSm4d(t,siglist,a2,d))
}
ellSm4dv3<-function(t,coord,siginv,ind){
  return(ellSm(c(pmax(t[ind[1]],t[ind[2]]),t[ind[3]]),coord[,2]-coord[,1],siginv))
}
ellSm4dv4<-function(t,coord,siginv){ 
  return(4*t[1]*t[2]*ellSm(t,coord[,2]-coord[,1],siginv))
}

ISmD2<-function(a){ 
  return(pnorm(a/2) + exp(a*a)*(1/3)*pnorm(-(3/2)*a))
}
I2SmD2<-function(x,a){ 
  lx <- log(x) 
  return(0.5*pnorm(a/2 - lx/a) + x*pnorm(a/2 + lx/a) + 0.5*x*x*exp(a*a)*pnorm(-(3/2)*a - lx/a))
} 
I3SmD2<-function(x,a){
  res<-.C("I3smith2",as.double(x),as.double(a),result=double(1),PACKAGE="spatialTailDep")$result
  return(res)
}

I1Sm<-function(u,a,siglist,a2,d,ind){ return(I3SmD2(u[ind],a)*ellSm4d(u,siglist,a2,d))} 
I1v2Sm<-function(u,a,siginv,coordt,ind,ind2){ return(I3SmD2(u[ind],a)*ellSm4dv3(u,coordt,siginv,ind2))} 
I1v3Sm<-function(u,a,siglist,a2,d){ return(I3SmD2(u[1],a)*ellSm4d(c(u[1],u[3],u[2]),siglist,a2,d))}
I1v4Sm<-function(u,a,siglist,a2,d){return(I3SmD2(u[1],a)*ellSm4d(c(u[2],u[1],u[3]),siglist,a2,d))} 
I3Sm<-function(u,coord,siginv,coordt){ 
  loc1<-coord[,2]-coord[,1]
  anr1<-sqrt((t(loc1)%*%siginv%*%loc1))
  loc2<-coord[,4]-coord[,3]
  anr2<-sqrt((t(loc2)%*%siginv%*%loc2))
  return(I3SmD2(u[1],anr1)*I3SmD2(u[2],anr2)*ellSm(u,coordt[,2]-coordt[,1],siginv))
}

## case 1: s,t,u,v all different
sigmaintSmD2case1<-function(coord,siginv,Tol){ 
  loc1<-coord[,2]-coord[,1]
  anr1<-sqrt((t(loc1)%*%siginv%*%loc1))
  loc2<-coord[,4]-coord[,3]
  anr2<-sqrt((t(loc2)%*%siginv%*%loc2))
  ########### T1
  d<-4
  siglist<-a2<-vector('list',length=d)
  for(i in 1:d){
    siglist[[i]]<-t(coord[,i]%*%t(rep(1,d-1)) - coord[,-i])%*%siginv%*%(coord[,i]%*%t(rep(1,d-1)) - coord[,-i])
    a2[[i]]<-diag(0.5*t(coord[,i]-coord[,-i])%*%siginv%*%(coord[,i]-coord[,-i]))
  }
  T1<-(ISmD2(anr1) + ISmD2(anr2) - adaptIntegrate(ellSm4d,lowerLimit=c(0,0,0,0),upperLimit=c(1,1,1,1),siglist=siglist,a2=a2,d=4,tol=Tol)$integral)
  ############# T2
  d<-3
  siglistpt1<-siglistpt2<-a2pt1<-a2pt2<-vector('list',length=3)
  coordt1<-coord[,-4]
  coordt2<-coord[,-3]
  for(i in 1:d){
    siglistpt1[[i]]<-t(coordt1[,i]%*%t(rep(1,d-1)) - coordt1[,-i])%*%siginv%*%(coordt1[,i]%*%t(rep(1,d-1)) - coordt1[,-i])
    a2pt1[[i]]<-diag(0.5*t(coordt1[,i]-coordt1[,-i])%*%siginv%*%(coordt1[,i]-coordt1[,-i]))
    siglistpt2[[i]]<-t(coordt2[,i]%*%t(rep(1,d-1)) - coordt2[,-i])%*%siginv%*%(coordt2[,i]%*%t(rep(1,d-1)) - coordt2[,-i])
    a2pt2[[i]]<-diag(0.5*t(coordt2[,i]-coordt2[,-i])%*%siginv%*%(coordt2[,i]-coordt2[,-i]))
  }
  T2<-(2*ISmD2(anr1)*(I2SmD2(1,anr2) -0.5) + 2*I2SmD2(1,anr2) - 2*ISmD2(anr2) - 
         adaptIntegrate(I1Sm,lowerLimit=c(0,0,0),upperLimit=c(1,1,1),a=anr2,siglist=siglistpt1,a2=a2pt1,d=3,ind=3,tol=Tol)$integral - 
         adaptIntegrate(I1Sm,lowerLimit=c(0,0,0),upperLimit=c(1,1,1),a=anr2,siglist=siglistpt2,a2=a2pt2,d=3,ind=3,tol=Tol)$integral)
  ############# T3
  d<-3
  siglistpt1<-siglistpt2<-a2pt1<-a2pt2<-vector('list',length=3)
  coordt1<-coord[,-2]
  coordt2<-coord[,-1]
  for(i in 1:d){
    siglistpt1[[i]]<-t(coordt1[,i]%*%t(rep(1,d-1)) - coordt1[,-i])%*%siginv%*%(coordt1[,i]%*%t(rep(1,d-1)) - coordt1[,-i])
    a2pt1[[i]]<-diag(0.5*t(coordt1[,i]-coordt1[,-i])%*%siginv%*%(coordt1[,i]-coordt1[,-i]))
    siglistpt2[[i]]<-t(coordt2[,i]%*%t(rep(1,d-1)) - coordt2[,-i])%*%siginv%*%(coordt2[,i]%*%t(rep(1,d-1)) - coordt2[,-i])
    a2pt2[[i]]<-diag(0.5*t(coordt2[,i]-coordt2[,-i])%*%siginv%*%(coordt2[,i]-coordt2[,-i]))
  }
  T3<-(2*ISmD2(anr2)*(I2SmD2(1,anr1) -0.5) + 2*I2SmD2(1,anr1) - 2*ISmD2(anr1) - 
         adaptIntegrate(I1Sm,lowerLimit=c(0,0,0),upperLimit=c(1,1,1),a=anr1,siglist=siglistpt1,a2=a2pt1,d=3,ind=1,tol=Tol)$integral -
         adaptIntegrate(I1Sm,lowerLimit=c(0,0,0),upperLimit=c(1,1,1),a=anr1,siglist=siglistpt2,a2=a2pt2,d=3,ind=1,tol=Tol)$integral)
  ################ T4
  T4<-((ISmD2(anr1) - I2SmD2(1,anr1))*(1 - 2*I2SmD2(1,anr2)) + (ISmD2(anr2) - I2SmD2(1,anr2))*(1 - 2*I2SmD2(1,anr1)) -
         adaptIntegrate(I3Sm,lowerLimit=c(0,0),upperLimit=c(1,1),coord=coord,siginv=siginv,coordt=coord[,c(1,3)])$integral - 
         adaptIntegrate(I3Sm,lowerLimit=c(0,0),upperLimit=c(1,1),coord=coord,siginv=siginv,coordt=coord[,c(2,4)])$integral)
  ################# T5
  T5<-((ISmD2(anr1) - I2SmD2(1,anr1))*(1 - 2*I2SmD2(1,anr2)) + (ISmD2(anr2) - I2SmD2(1,anr2))*(1 - 2*I2SmD2(1,anr1)) -
         adaptIntegrate(I3Sm,lowerLimit=c(0,0),upperLimit=c(1,1),coord=coord,siginv=siginv,coordt=coord[,c(1,4)])$integral - 
         adaptIntegrate(I3Sm,lowerLimit=c(0,0),upperLimit=c(1,1),coord=coord,siginv=siginv,coordt=coord[,c(2,3)])$integral)
  return(T1 - T2 - T3 + T4 + T5)
}

## case 2: t = u
sigmaintSmD2case2<-function(coord,siginv,Tol){ 
  loc1<-coord[,2]-coord[,1]
  anr1<-sqrt((t(loc1)%*%siginv%*%loc1))
  loc2<-coord[,4]-coord[,3]
  anr2<-sqrt((t(loc2)%*%siginv%*%loc2))
  ########### T1
  d<-3
  siglist<-a2<-vector('list',length=d)
  coordt<-coord[,-3]
  for(i in 1:d){
    siglist[[i]]<-t(coordt[,i]%*%t(rep(1,d-1)) - coordt[,-i])%*%siginv%*%(coordt[,i]%*%t(rep(1,d-1)) - coordt[,-i])
    a2[[i]]<-diag(0.5*t(coordt[,i]-coordt[,-i])%*%siginv%*%(coordt[,i]-coordt[,-i]))
  }
  T1<-(ISmD2(anr1) + ISmD2(anr2) - 
         adaptIntegrate(ellSm4dv2,lowerLimit=c(0,0,0),upperLimit=c(1,1,1),siglist=siglist,a2=a2,d=3,ind=2,tol=Tol)$integral)
  ############# T2
  T2<-(2*ISmD2(anr1)*(I2SmD2(1,anr2) -0.5) + 2*I2SmD2(1,anr2) - 2*ISmD2(anr2) - 
         adaptIntegrate(I1v2Sm,lowerLimit=c(0,0,0),upperLimit=c(1,1,1),a=anr2,siginv=siginv,
                        coordt=coord[,1:2],ind=3,ind2=c(2,3,1),tol=Tol)$integral - 
         adaptIntegrate(I1Sm,lowerLimit=c(0,0,0),upperLimit=c(1,1,1),a=anr2,siglist=siglist,a2=a2,d=3,ind=3,tol=Tol)$integral)
  ############# T3
  T3<-(2*ISmD2(anr2)*(I2SmD2(1,anr1) -0.5) + 2*I2SmD2(1,anr1) - 2*ISmD2(anr1) - 
         adaptIntegrate(I1Sm,lowerLimit=c(0,0,0),upperLimit=c(1,1,1),a=anr1,siglist=siglist,a2=a2,d=3,ind=1,tol=Tol)$integral - 
         adaptIntegrate(I1v2Sm,lowerLimit=c(0,0,0),upperLimit=c(1,1,1),a=anr1,siginv=siginv,
                        coordt=coord[,3:4],ind=1,ind2=c(1,2,3),tol=Tol)$integral)    
  ################ T4
  T4<-((ISmD2(anr1) - I2SmD2(1,anr1))*(1 - 2*I2SmD2(1,anr2)) + (ISmD2(anr2) - I2SmD2(1,anr2))*(1 - 2*I2SmD2(1,anr1)) -
         adaptIntegrate(I3Sm,lowerLimit=c(0,0),upperLimit=c(1,1),coord=coord,siginv=siginv,coordt=coord[,c(1,3)])$integral - 
         adaptIntegrate(I3Sm,lowerLimit=c(0,0),upperLimit=c(1,1),coord=coord,siginv=siginv,coordt=coord[,c(2,4)])$integral)
  ################# T5
  T5<-((ISmD2(anr1) - I2SmD2(1,anr1))*(1 - 2*I2SmD2(1,anr2)) + (ISmD2(anr2) - I2SmD2(1,anr2))*(1 - 2*I2SmD2(1,anr1)) -
         adaptIntegrate(I3Sm,lowerLimit=c(0,0),upperLimit=c(1,1),coord=coord,siginv=siginv,coordt=coord[,c(1,4)])$integral - 
         adaptIntegrate(I3Sm,lowerLimit=c(0,0),upperLimit=c(1,1),coord=coord,siginv=siginv,coordt=coord[,c(2,3)])$integral)
  return(T1 - T2 - T3 + T4 + T5)
}

## case 3: s = u, t neq v
sigmaintSmD2case3<-function(coord,siginv,Tol){ 
  loc1<-coord[,2]-coord[,1]
  anr1<-sqrt((t(loc1)%*%siginv%*%loc1))
  loc2<-coord[,4]-coord[,3]
  anr2<-sqrt((t(loc2)%*%siginv%*%loc2))
  ########### T1
  d<-3
  siglist<-a2<-vector('list',length=d)
  coordt<-coord[,-3]
  for(i in 1:d){
    siglist[[i]]<-t(coordt[,i]%*%t(rep(1,d-1)) - coordt[,-i])%*%siginv%*%(coordt[,i]%*%t(rep(1,d-1)) - coordt[,-i])
    a2[[i]]<-diag(0.5*t(coordt[,i]-coordt[,-i])%*%siginv%*%(coordt[,i]-coordt[,-i]))
  }
  T1<-(ISmD2(anr1) + ISmD2(anr2) - 
         adaptIntegrate(ellSm4dv2,lowerLimit=c(0,0,0),upperLimit=c(1,1,1),siglist=siglist,a2=a2,d=3,ind=1,tol=Tol)$integral)
  ############# T2
  T2<-(2*ISmD2(anr1)*(I2SmD2(1,anr2) -0.5) + 2*I2SmD2(1,anr2) - 2*ISmD2(anr2) - 
         adaptIntegrate(I1v2Sm,lowerLimit=c(0,0,0),upperLimit=c(1,1,1),a=anr2,siginv=siginv,
                        coordt=coord[,1:2],ind=3,ind2=c(1,3,2),tol=Tol)$integral - 
         adaptIntegrate(I1Sm,lowerLimit=c(0,0,0),upperLimit=c(1,1,1),a=anr2,siglist=siglist,a2=a2,d=3,ind=3,tol=Tol)$integral)
  ############# T3
  T3<-(2*ISmD2(anr2)*(I2SmD2(1,anr1) -0.5) + 2*I2SmD2(1,anr1) - 2*ISmD2(anr1) - 
         adaptIntegrate(I1v2Sm,lowerLimit=c(0,0,0),upperLimit=c(1,1,1),a=anr1,siginv=siginv,
                        coordt=coord[,c(1,4)],ind=1,ind2=c(1,2,3),tol=Tol)$integral - 
         adaptIntegrate(I1v4Sm,lowerLimit=c(0,0,0),upperLimit=c(1,1,1), a=anr1,siglist=siglist,a2=a2,d=3,tol=Tol)$integral)   
  ################ T4
  T4<-((ISmD2(anr1) - I2SmD2(1,anr1))*(1 - 2*I2SmD2(1,anr2)) + (ISmD2(anr2) - I2SmD2(1,anr2))*(1 - 2*I2SmD2(1,anr1)) -
         adaptIntegrate(I3Sm,lowerLimit=c(0,0),upperLimit=c(1,1),coord=coord,siginv=siginv,coordt=coord[,c(1,3)])$integral - 
         adaptIntegrate(I3Sm,lowerLimit=c(0,0),upperLimit=c(1,1),coord=coord,siginv=siginv,coordt=coord[,c(2,4)])$integral)
  ################# T5
  T5<-((ISmD2(anr1) - I2SmD2(1,anr1))*(1 - 2*I2SmD2(1,anr2)) + (ISmD2(anr2) - I2SmD2(1,anr2))*(1 - 2*I2SmD2(1,anr1)) -
         adaptIntegrate(I3Sm,lowerLimit=c(0,0),upperLimit=c(1,1),coord=coord,siginv=siginv,coordt=coord[,c(1,4)])$integral - 
         adaptIntegrate(I3Sm,lowerLimit=c(0,0),upperLimit=c(1,1),coord=coord,siginv=siginv,coordt=coord[,c(2,3)])$integral)
  return(T1 - T2 - T3 + T4 + T5)
}

## case 4: s neq u, t = v
sigmaintSmD2case4<-function(coord,siginv,Tol){ 
  loc1<-coord[,2]-coord[,1]
  anr1<-sqrt((t(loc1)%*%siginv%*%loc1))
  loc2<-coord[,4]-coord[,3]
  anr2<-sqrt((t(loc2)%*%siginv%*%loc2))
  ########### T1
  d<-3
  siglist<-a2<-vector('list',length=d)
  coordt<-coord[,-4]
  for(i in 1:d){
    siglist[[i]]<-t(coordt[,i]%*%t(rep(1,d-1)) - coordt[,-i])%*%siginv%*%(coordt[,i]%*%t(rep(1,d-1)) - coordt[,-i])
    a2[[i]]<-diag(0.5*t(coordt[,i]-coordt[,-i])%*%siginv%*%(coordt[,i]-coordt[,-i]))
  }
  T1<-(ISmD2(anr1) + ISmD2(anr2) - 
         adaptIntegrate(ellSm4dv2,lowerLimit=c(0,0,0),upperLimit=c(1,1,1),siglist=siglist,a2=a2,d=3,ind=2,tol=Tol)$integral)
  ############# T2
  T2<-(2*ISmD2(anr1)*(I2SmD2(1,anr2) -0.5) + 2*I2SmD2(1,anr2) - 2*ISmD2(anr2) - 
         adaptIntegrate(I1Sm,lowerLimit=c(0,0,0),upperLimit=c(1,1,1),a=anr2,siglist=siglist,a2=a2,d=3,ind=3,tol=Tol)$integral -
         adaptIntegrate(I1v2Sm,lowerLimit=c(0,0,0),upperLimit=c(1,1,1),a=anr2,siginv=siginv,
                        coordt=coord[,c(1,2)],ind=3,ind2=c(2,3,1),tol=Tol)$integral) 
  ############# T3
  T3<-(2*ISmD2(anr2)*(I2SmD2(1,anr1) -0.5) + 2*I2SmD2(1,anr1) - 2*ISmD2(anr1) - 
         adaptIntegrate(I1v3Sm,lowerLimit=c(0,0,0),upperLimit=c(1,1,1),a=anr1,siglist=siglist,a2=a2,d=3,tol=Tol)$integral - 
         adaptIntegrate(I1v2Sm,lowerLimit=c(0,0,0),upperLimit=c(1,1,1),a=anr1,siginv=siginv,
                        coordt=coord[,c(2,3)],ind=1,ind2=c(1,3,2),tol=Tol)$integral)    
  ################ T4
  T4<-((ISmD2(anr1) - I2SmD2(1,anr1))*(1 - 2*I2SmD2(1,anr2)) + (ISmD2(anr2) - I2SmD2(1,anr2))*(1 - 2*I2SmD2(1,anr1)) -
         adaptIntegrate(I3Sm,lowerLimit=c(0,0),upperLimit=c(1,1),coord=coord,siginv=siginv,coordt=coord[,c(1,3)])$integral - 
         adaptIntegrate(I3Sm,lowerLimit=c(0,0),upperLimit=c(1,1),coord=coord,siginv=siginv,coordt=coord[,c(2,4)])$integral)
  ################# T5
  T5<-((ISmD2(anr1) - I2SmD2(1,anr1))*(1 - 2*I2SmD2(1,anr2)) + (ISmD2(anr2) - I2SmD2(1,anr2))*(1 - 2*I2SmD2(1,anr1)) -
         adaptIntegrate(I3Sm,lowerLimit=c(0,0),upperLimit=c(1,1),coord=coord,siginv=siginv,coordt=coord[,c(1,4)])$integral - 
         adaptIntegrate(I3Sm,lowerLimit=c(0,0),upperLimit=c(1,1),coord=coord,siginv=siginv,coordt=coord[,c(2,3)])$integral)
  return(T1 - T2 - T3 + T4 + T5)
}

## case 5: s = u, t = v
sigmaintSmD2case5<-function(coord,siginv,Tol){ 
  loc<-coord[,2]-coord[,1]
  a<-sqrt((t(loc)%*%siginv%*%loc))
  ########### T1
  coordt<-coord[,1:2]
  T1<-(2*ISmD2(a)  - adaptIntegrate(ellSm4dv4,lowerLimit=c(0,0),upperLimit=c(1,1),coord=coordt,siginv=siginv,tol=Tol)$integral)
  ############# T2
  T2<-(2*ISmD2(a)*(I2SmD2(1,a) -0.5) + 2*I2SmD2(1,a) - 2*ISmD2(a) - 
         adaptIntegrate(I1v2Sm,lowerLimit=c(0,0,0),upperLimit=c(1,1,1),a=a,siginv=siginv,
                        coordt=coordt,ind=3,ind2=c(1,3,2),tol=Tol)$integral - 
         adaptIntegrate(I1v2Sm,lowerLimit=c(0,0,0),upperLimit=c(1,1,1),a=a,siginv=siginv,
                        coordt=coordt,ind=3,ind2=c(2,3,1),tol=Tol)$integral)
  ############# T3
  T3<-(2*ISmD2(a)*(I2SmD2(1,a) -0.5) + 2*I2SmD2(1,a) - 2*ISmD2(a) - 
         adaptIntegrate(I1v2Sm,lowerLimit=c(0,0,0),upperLimit=c(1,1,1),a=a,siginv=siginv,
                        coordt=coordt,ind=1,ind2=c(1,2,3),tol=Tol)$integral - 
         adaptIntegrate(I1v2Sm,lowerLimit=c(0,0,0),upperLimit=c(1,1,1),a=a,siginv=siginv,
                        coordt=coordt,ind=1,ind2=c(1,3,2),tol=Tol)$integral)    
  ################ T4
  T4<-((ISmD2(a) - I2SmD2(1,a))*(1 - 2*I2SmD2(1,a)) + (ISmD2(a) - I2SmD2(1,a))*(1 - 2*I2SmD2(1,a)) -
         adaptIntegrate(I3Sm,lowerLimit=c(0,0),upperLimit=c(1,1),coord=coord,siginv=siginv,coordt=coord[,c(1,3)])$integral - 
         adaptIntegrate(I3Sm,lowerLimit=c(0,0),upperLimit=c(1,1),coord=coord,siginv=siginv,coordt=coord[,c(2,4)])$integral)
  ################# T5
  T5<-((ISmD2(a) - I2SmD2(1,a))*(1 - 2*I2SmD2(1,a)) + (ISmD2(a) - I2SmD2(1,a))*(1 - 2*I2SmD2(1,a)) -
         adaptIntegrate(I3Sm,lowerLimit=c(0,0),upperLimit=c(1,1),coord=coord,siginv=siginv,coordt=coord[,c(1,4)])$integral - 
         adaptIntegrate(I3Sm,lowerLimit=c(0,0),upperLimit=c(1,1),coord=coord,siginv=siginv,coordt=coord[,c(2,3)])$integral)
  return(T1 - T2 - T3 + T4 + T5)
}

psidotSmD2<-function(sigmam,pairs){
  siginv<-solve(sigmam)
  npairs<-nrow(pairs)
  loc<-pairs[,3:4]-pairs[,1:2]
  psi<-matrix(,nrow=npairs,ncol=3)
  noem<-(sigmam[1,1]*sigmam[2,2] - sigmam[1,2]*sigmam[2,1])
  noem2 <- noem * noem 
  deriv<-vector(length=3)
  for(i in 1:npairs){
    a2 <- t(loc[i,]) %*% siginv %*% loc[i,]
    a <- sqrt(a2)
    tel<-(sigmam[2,2]*((loc[i,1])^2) + sigmam[1,1]*((loc[i,2])^2) - 
            2*sigmam[1,2]*(loc[i,1])*(loc[i,2]))
    deriv[1] <- (((loc[i,2])^2)*noem - sigmam[2,2]*tel)
    deriv[2] <- (((loc[i,1])^2)*noem - sigmam[1,1]*tel)
    deriv[3] <- (2*tel*sigmam[1,2] - 2*noem*loc[i,1]*loc[i,2])
    deriv <- deriv / (2 * a * noem2)
    psi[i,]<-deriv*((1/2)*dnorm(a/2) - (1/2)*exp(a2)*dnorm(-(3/2)*a) + (2/3)*a*exp(a2)*pnorm(-(3/2)*a))
  }
  return(psi)
}


AsymVarSm<-function(pairs, sigmainv, Tol){
  npairs<-nrow(pairs)
  Mat<-matrix(,nrow=npairs,ncol=npairs)
  for(i in 1:npairs){
    for(j in i:npairs){
      temp<-matrix(c(pairs[i,],pairs[j,]),nrow=2)
      if(all(temp[,1]==temp[,4]) && (! all(temp[,2]==temp[,3]))){
        Mat[i,j]<-sigmaintSmD2case2(cbind(temp[,3:4],temp[,1:2]),sigmainv,Tol=Tol)
      }else if(all(temp[,1]==temp[,3]) && all(temp[,2]==temp[,4])){
        Mat[i,j]<-sigmaintSmD2case5(temp,sigmainv,Tol=Tol)
      }else if(all(temp[,2]==temp[,4]) && (! all(temp[,1]==temp[,3]))){
        Mat[i,j]<-sigmaintSmD2case4(temp,sigmainv,Tol=Tol)
      }else if(all(temp[,1]==temp[,3]) && (! all(temp[,2]==temp[,4]))){
        Mat[i,j]<-sigmaintSmD2case3(temp,sigmainv,Tol=Tol)
      }else if(all(temp[,2]==temp[,3]) && (! all(temp[,1]==temp[,4]))){
        Mat[i,j]<-sigmaintSmD2case2(temp,sigmainv,Tol=Tol)
      }else if((! all(temp[,1]==temp[,3])) && (! all(temp[,2]==temp[,4]))){
        Mat[i,j]<-sigmaintSmD2case1(temp,sigmainv,Tol=Tol)
      }
      Mat[j,i]<-Mat[i,j]
    }
  }
  return(Mat)
}

##### change to l-bfgs and control argument
MestimatorSmith <- function(x, locations, pairIndices, k,Tol, startingValue, Omega, iterate, covMat ) {
  n <- dim(x)[1]
  if (k < 1 || k > (n-1)) {
    warning("k must be between 1 and n")
  }
  npairs <- nrow(pairIndices)
  pairs <- pairCoordinates(locations, pairIndices)
  ranks <- apply(x, 2, rank) 
  totL <- numeric(npairs) 
  for (l in 1:npairs) {
    pairsOfRanks <- cbind(ranks[, pairIndices[l, 1]],ranks[, pairIndices[l, 2]])
    totL[l] <- tailIntEmp(pairsOfRanks, n, k)
  }
  strt<-c(startingValue[1,1],startingValue[2,2],startingValue[1,2])
  sigma <- optim(par = strt, fn = tominimizeSmith, totlist = totL,pairs = pairs,w = Omega)$par
  sigmam <- sigmamPilot <- rbind(c(sigma[1], sigma[3]),c(sigma[3], sigma[2]))
  siginv <- solve(sigmam, nrow=2)
   if (iterate) {
    Omega <- solve(AsymVarSm(pairs, siginv, Tol))
    sigma <- optim(par = c(sigmamPilot[1,1],sigmamPilot[2,2],sigmamPilot[1,2]),fn = tominimizeSmith,totlist = totL,pairs = pairs,w = Omega)$par    
    sigmam <- rbind(c(sigma[1], sigma[3]), c(sigma[3], sigma[2]))
    siginv <- solve(sigmam, nrow=2)
  }
  covMatrix<-NULL
  if(covMat){
    psidot <- psidotSmD2(sigmam, pairs)
    temp <- solve(t(psidot) %*% Omega %*% psidot) %*% t(psidot)
    GammaMat <- AsymVarSm(pairs, siginv, Tol)
    covMatrix <- (temp %*% Omega %*% GammaMat %*% Omega %*% t(temp)) / k
  }
  return(list(theta = sigmam,theta_pilot = sigmamPilot,covMatrix = covMatrix,Omega = Omega))
}
