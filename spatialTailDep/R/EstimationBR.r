#' @include gridUtilities.r
#' @include tailIntEmp.r
NULL

parsTauToBR<-function(tau){
  wortel<-sqrt(1 + (4*(tau[3]^2))/((tau[2]-tau[1])^2))
  rat<-(tau[2]+tau[1])/(tau[1]-tau[2])
  beta<-0.5*atan((2*tau[3])/(tau[2]-tau[1]))
  if(beta >=0){
    cc<-sqrt((rat-wortel)/(rat+wortel))
    rho<-sqrt((cc^2 - 1)/(wortel*(tau[2]-tau[1])))
  } else{
    beta<-0.5*(atan((2*tau[3])/(tau[2]-tau[1]))+pi)
    cc<-1/sqrt((rat-wortel)/(rat+wortel))
    rho<-sqrt((cc^2 - 1)/(wortel*(tau[1]-tau[2])))
  }  
  return(c(rho,beta,cc))
}

parsBRtoTau<-function(pars){
  tau11<-(1/(2*(pars[1]^2)))*(1 + pars[3]^2 - (pars[3]^2 - 1)*cos(2*pars[2]))
  tau22<-(1/(2*(pars[1]^2)))*(1 + pars[3]^2 + (pars[3]^2 - 1)*cos(2*pars[2]))
  tau12<-(((pars[3]^2 - 1)/(2*(pars[1]^2)))*sin(2*pars[2]))
  return(c(tau11,tau22,tau12))
}

tailBRInt <- function(loc, tau, alpha){
  locdiff<-loc[3:4]-loc[1:2]
  a <- sqrt(2)*((t(locdiff) %*% tau %*% locdiff)^(alpha/4))
  logres<-(a^2 + pnorm((-3*a)/2,log.p=TRUE) - log(3))
  return(as.vector(pnorm(a/2) + exp(logres)))
}

tominimizeBR<-function(pars,totlist,pairs,w){
  tau <- rbind(c(pars[2], pars[4]), c(pars[4], pars[3]))
  res<-apply(pairs,1,tailBRInt,tau=tau,alpha=pars[1]) - totlist
  return(t(res) %*% w %*% res)
}


ellBR<-function(t,loc,tau,alpha){ 
  if((loc[1]==0) && (loc[2]==0)){
    return(pmax(t[1],t[2]))
  } else{
    locdiff<-loc[3:4]-loc[1:2]
    a <- sqrt(2)*((t(loc) %*% tau %*% loc)^(alpha/4))
    res<-.C("ellsmith2",as.double(t),as.double(a),result=double(1),PACKAGE="spatialTailDep")$result
    return(res)
  }
}

ellBR4d<-function(t,rlist,gamlist,d){ 
  result<-vector(length=d)
  for(i in 1:d){
    const<-(gamlist[[i]]/2 + log(t[i]/t[-i])/gamlist[[i]])
    result[i]<-t[i]*pmvnorm(lower=rep(-Inf,d-1),upper=const,corr=rlist[[i]],algorithm=TVPACK)[1]
  }
  return(sum(result))
}
gfunc<-function(tau,alpha,loc){ 
  a2 <- ((t(loc) %*% tau %*% loc)^(alpha/2))
  return(a2)
}

corrfunc<-function(tau,alpha,loc,d,i){ 
  result<-matrix(1,nrow=(d-1),ncol=(d-1))
  for(j in 1:(d-2)){
    for(k in (j+1):(d-1)){
      lo<-loc[,-i]
      result[j,k]<-result[k,j]<-((gfunc(tau,alpha,loc[,i]-lo[,j]) + gfunc(tau,alpha,loc[,i]-lo[,k]) -
                                    gfunc(tau,alpha,lo[,j]-lo[,k]))/(2*sqrt(gfunc(tau,alpha,loc[,i]-lo[,j])*gfunc(tau,alpha,loc[,i]-lo[,k]))))
    }
  }
  return(result)
}

tailBR <- function(t, loc, tau, alpha, d){ 
  if(d==2){
    return(ellBR(t,loc[,2]-loc[,1],tau,alpha))
  } else{
    rlist<-gamlist<-vector('list',length=d)
    for(i in 1:d){
      rlist[[i]]<-corrfunc(tau,alpha,loc,d,i)
      gamlist[[i]]<-sapply(c(1:d)[-i], function(j) sqrt(2*gfunc(tau,alpha,loc[,i]-loc[,j])))
    }
    return(ellBR4d(t,rlist,gamlist,d))
  }
}

ellBR4dv2<-function(t,rlist,gamlist,d,ind){ 
  return(2*t[ind]*ellBR4d(t,rlist,gamlist,d))
}
ellBR4dv3<-function(t,coord,tau,alpha,ind){
  return(ellBR(c(pmax(t[ind[1]],t[ind[2]]),t[ind[3]]),coord[,2]-coord[,1],tau,alpha))
}
ellBR4dv4<-function(t,coord,tau,alpha){ 
  return(4*t[1]*t[2]*ellBR(t,coord[,2]-coord[,1],tau,alpha))
}

IBR<-function(a){ #  This function can be implemented in C as well.
  return(pnorm(a/2) + exp(a*a)*(1/3)*pnorm(-(3/2)*a))
}
I2BR<-function(x,a){  #  This function can be implemented in C as well.
  return(0.5*pnorm(a/2 - log(x)/a) + x*pnorm(a/2 + log(x)/a) + 0.5*x*x*exp(a*a)*pnorm(-(3/2)*a - log(x)/a))
} 
I3BR<-function(x,a){
  res<-.C("I3smith2",as.double(x),as.double(a),result=double(1),PACKAGE="spatialTailDep")$result
  return(res)
}

I1<-function(u,a,rlist,gamlist,d,ind){ return(I3BR(u[ind],a)*ellBR4d(u,rlist,gamlist,d))} 
I1v2<-function(u,a,tau,alpha,coordt,ind,ind2){ return(I3BR(u[ind],a)*ellBR4dv3(u,coordt,tau,alpha,ind2))} 
I1v3<-function(u,a,rlist,gamlist,d){return(I3BR(u[1],a)*ellBR4d(c(u[1],u[3],u[2]),rlist,gamlist,d))} 
I1v4<-function(u,a,rlist,gamlist,d){return(I3BR(u[1],a)*ellBR4d(c(u[2],u[1],u[3]),rlist,gamlist,d))} 
I3<-function(u,coord,tau,alpha,coordt){ 
  loc1<-coord[,2]-coord[,1]
  loc2<-coord[,4]-coord[,3]
  anr1 <- sqrt(2)*((t(loc1) %*% tau %*% loc1)^(alpha/4))
  anr2 <- sqrt(2)*((t(loc2) %*% tau %*% loc2)^(alpha/4))
  return(I3BR(u[1],anr1)*I3BR(u[2],anr2)*ellBR(u,coordt[,2]-coordt[,1],tau,alpha))}


####### case 1: locations s,t,u,v all different
intBRcase1<-function(coord,tau,alpha,Tol){ #coord=(s,t,u,v) where s,t,u,v, are 2-dim columnvectors
  loc1<-coord[,2]-coord[,1]
  loc2<-coord[,4]-coord[,3]
  anr1 <- sqrt(2)*((t(loc1) %*% tau %*% loc1)^(alpha/4))
  anr2 <- sqrt(2)*((t(loc2) %*% tau %*% loc2)^(alpha/4))
  ########### T1
  d<-4
  rlist<-gamlist<-vector('list',length=d)
  for(i in 1:d){
    rlist[[i]]<-corrfunc(tau,alpha,coord,d,i)
    gamlist[[i]]<-sapply(c(1:d)[-i], function(j) sqrt(2*gfunc(tau,alpha,coord[,i]-coord[,j])))
  }
  T1<-(IBR(anr1) + IBR(anr2) - adaptIntegrate(ellBR4d,lowerLimit=c(0,0,0,0),upperLimit=c(1,1,1,1),rlist=rlist,gamlist=gamlist,d=4,tol=Tol)$integral)
  ############# T2
  d<-3
  rlistpt1<-rlistpt2<-gamlistpt1<-gamlistpt2<-vector('list',length=d)
  coordt1<-coord[,-4]
  coordt2<-coord[,-3]
  for(i in 1:d){
    rlistpt1[[i]]<-corrfunc(tau,alpha,coordt1,d,i)
    gamlistpt1[[i]]<-sapply(c(1:d)[-i], function(j) sqrt(2*gfunc(tau,alpha,coordt1[,i]-coordt1[,j])))
    rlistpt2[[i]]<-corrfunc(tau,alpha,coordt2,d,i)
    gamlistpt2[[i]]<-sapply(c(1:d)[-i], function(j) sqrt(2*gfunc(tau,alpha,coordt2[,i]-coordt2[,j])))
  }
  T2<-(2*IBR(anr1)*(I2BR(1,anr2) -0.5) + 2*I2BR(1,anr2) - 2*IBR(anr2) - 
         adaptIntegrate(I1,lowerLimit=c(0,0,0),upperLimit=c(1,1,1),a=anr2,rlist=rlistpt1,gamlist=gamlistpt1,d=3,ind=3,tol=Tol)$integral - 
         adaptIntegrate(I1,lowerLimit=c(0,0,0),upperLimit=c(1,1,1),a=anr2,rlist=rlistpt2,gamlist=gamlistpt2,d=3,ind=3,tol=Tol)$integral)
  ############# T3
  d<-3
  rlistpt1<-rlistpt2<-gamlistpt1<-gamlistpt2<-vector('list',length=3)
  coordt1<-coord[,-2]
  coordt2<-coord[,-1]
  for(i in 1:d){
    rlistpt1[[i]]<-corrfunc(tau,alpha,coordt1,d,i)
    gamlistpt1[[i]]<-sapply(c(1:d)[-i], function(j) sqrt(2*gfunc(tau,alpha,coordt1[,i]-coordt1[,j])))
    rlistpt2[[i]]<-corrfunc(tau,alpha,coordt2,d,i)
    gamlistpt2[[i]]<-sapply(c(1:d)[-i], function(j) sqrt(2*gfunc(tau,alpha,coordt2[,i]-coordt2[,j])))
  }
  T3<-(2*IBR(anr2)*(I2BR(1,anr1) -0.5) + 2*I2BR(1,anr1) - 2*IBR(anr1) - 
         adaptIntegrate(I1,lowerLimit=c(0,0,0),upperLimit=c(1,1,1), a=anr1,rlist=rlistpt1,gamlist=gamlistpt1,d=3,ind=1,tol=Tol)$integral -
         adaptIntegrate(I1,lowerLimit=c(0,0,0),upperLimit=c(1,1,1), a=anr1,rlist=rlistpt2,gamlist=gamlistpt2,d=3,ind=1,tol=Tol)$integral)
  ################ T4
  T4<-((IBR(anr1) - I2BR(1,anr1))*(1 - 2*I2BR(1,anr2)) + (IBR(anr2) -  I2BR(1,anr2))*(1 - 2*I2BR(1,anr1)) -
         adaptIntegrate(I3,lowerLimit=c(0,0),upperLimit=c(1,1),coord=coord,tau=tau,alpha=alpha,coordt=coord[,c(1,3)])$integral - 
         adaptIntegrate(I3,lowerLimit=c(0,0),upperLimit=c(1,1),coord=coord,tau=tau,alpha=alpha,coordt=coord[,c(2,4)])$integral)
  ################# T5
  T5<-((IBR(anr1) - I2BR(1,anr1))*(1 - 2*I2BR(1,anr2)) + (IBR(anr2) -  I2BR(1,anr2))*(1 - 2*I2BR(1,anr1)) -
         adaptIntegrate(I3,lowerLimit=c(0,0),upperLimit=c(1,1),coord=coord,tau=tau,alpha=alpha,coordt=coord[,c(1,4)])$integral - 
         adaptIntegrate(I3,lowerLimit=c(0,0),upperLimit=c(1,1),coord=coord,tau=tau,alpha=alpha,coordt=coord[,c(2,3)])$integral)
  return(T1 - T2 - T3 + T4 + T5)
}

## case 2: t = u
intBRcase2<-function(coord,tau,alpha,Tol){ #coord=(s,t,u,v) where s,t,u,v, are 2-dim columnvectors
  loc1<-coord[,2]-coord[,1]
  loc2<-coord[,4]-coord[,3]
  anr1 <- sqrt(2)*((t(loc1) %*% tau %*% loc1)^(alpha/4))
  anr2 <- sqrt(2)*((t(loc2) %*% tau %*% loc2)^(alpha/4))
  ########### T1
  d<-3
  rlist<-gamlist<-vector('list',length=d)
  coordt<-coord[,-3]
  for(i in 1:d){
    rlist[[i]]<-corrfunc(tau,alpha,coordt,d,i)
    gamlist[[i]]<-sapply(c(1:d)[-i], function(j) sqrt(2*gfunc(tau,alpha,coordt[,i]-coordt[,j])))
  }  
  T1<-(IBR(anr1) + IBR(anr2) - 
         adaptIntegrate(ellBR4dv2,lowerLimit=c(0,0,0),upperLimit=c(1,1,1),rlist=rlist,gamlist=gamlist,d=3,ind=2,tol=Tol)$integral)
  ############# T2
  T2<-(2*IBR(anr1)*(I2BR(1,anr2) -0.5) + 2*I2BR(1,anr2) - 2*IBR(anr2) - 
         adaptIntegrate(I1v2,lowerLimit=c(0,0,0),upperLimit=c(1,1,1),a=anr2,tau=tau,alpha=alpha,coordt=coord[,1:2],ind=3,ind2=c(2,3,1),tol=Tol)$integral - 
         adaptIntegrate(I1,lowerLimit=c(0,0,0),upperLimit=c(1,1,1),a=anr2,rlist=rlist,gamlist=gamlist,d=3,ind=3,tol=Tol)$integral)
  ############# T3
  T3<-(2*IBR(anr2)*(I2BR(1,anr1) -0.5) + 2*I2BR(1,anr1) - 2*IBR(anr1) - 
         adaptIntegrate(I1,lowerLimit=c(0,0,0),upperLimit=c(1,1,1),a=anr1, rlist=rlist,gamlist=gamlist,d=3,ind=1,tol=Tol)$integral - 
         adaptIntegrate(I1v2,lowerLimit=c(0,0,0),upperLimit=c(1,1,1),a=anr1,tau=tau,alpha=alpha,coordt=coord[,3:4],ind=1,ind2=c(1,2,3),tol=Tol)$integral)    
  ################ T4
  T4<-((IBR(anr1) - I2BR(1,anr1))*(1 - 2*I2BR(1,anr2)) + (IBR(anr2) -  I2BR(1,anr2))*(1 - 2*I2BR(1,anr1)) -
         adaptIntegrate(I3,lowerLimit=c(0,0),upperLimit=c(1,1),coord=coord,tau=tau,alpha=alpha,coordt=coord[,c(1,3)])$integral - 
         adaptIntegrate(I3,lowerLimit=c(0,0),upperLimit=c(1,1),coord=coord,tau=tau,alpha=alpha,coordt=coord[,c(2,4)])$integral)
  ################# T5
  T5<-((IBR(anr1) - I2BR(1,anr1))*(1 - 2*I2BR(1,anr2)) + (IBR(anr2) -  I2BR(1,anr2))*(1 - 2*I2BR(1,anr1)) -
         adaptIntegrate(I3,lowerLimit=c(0,0),upperLimit=c(1,1),coord=coord,tau=tau,alpha=alpha,coordt=coord[,c(1,4)])$integral - 
         adaptIntegrate(I3,lowerLimit=c(0,0),upperLimit=c(1,1),coord=coord,tau=tau,alpha=alpha,coordt=coord[,c(2,3)])$integral)
  return(T1 - T2 - T3 + T4 + T5)
}

## case 3: s = u, t neq v
intBRcase3<-function(coord,tau,alpha,Tol){ #coord=(s,t,u,v) where s,t,u,v, are 2-dim columnvectors
  loc1<-coord[,2]-coord[,1]
  loc2<-coord[,4]-coord[,3]
  anr1 <- sqrt(2)*((t(loc1) %*% tau %*% loc1)^(alpha/4))
  anr2 <- sqrt(2)*((t(loc2) %*% tau %*% loc2)^(alpha/4))
  ########### T1
  d<-3
  rlist<-gamlist<-vector('list',length=d)
  coordt<-coord[,-3]
  for(i in 1:d){
    rlist[[i]]<-corrfunc(tau,alpha,coordt,d,i)
    gamlist[[i]]<-sapply(c(1:d)[-i], function(j) sqrt(2*gfunc(tau,alpha,coordt[,i]-coordt[,j])))
  }
  T1<-(IBR(anr1)+IBR(anr2)-adaptIntegrate(ellBR4dv2,lowerLimit=c(0,0,0),upperLimit=c(1,1,1),rlist=rlist,gamlist=gamlist,d=3,ind=1,tol=Tol)$integral)
  ############# T2
  T2<-(2*IBR(anr1)*(I2BR(1,anr2) - 0.5) + 2*I2BR(1,anr2) - 2*IBR(anr2) - 
         adaptIntegrate(I1v2,lowerLimit=c(0,0,0),upperLimit=c(1,1,1),a=anr2,tau=tau,alpha=alpha,coordt=coord[,1:2],ind=3,ind2=c(1,3,2),tol=Tol)$integral - 
         adaptIntegrate(I1,lowerLimit=c(0,0,0),upperLimit=c(1,1,1),a=anr2,rlist=rlist,gamlist=gamlist,d=3,ind=3,tol=Tol)$integral)
  ############# T3
  T3<-(2*IBR(anr2)*(I2BR(1,anr1) -0.5) + 2*I2BR(1,anr1) - 2*IBR(anr1) - 
         adaptIntegrate(I1v2,lowerLimit=c(0,0,0),upperLimit=c(1,1,1),a=anr1,tau=tau,alpha=alpha,coordt=coord[,c(1,4)],ind=1,ind2=c(1,2,3),tol=Tol)$integral - 
         adaptIntegrate(I1v4,lowerLimit=c(0,0,0),upperLimit=c(1,1,1),a=anr1,rlist=rlist,gamlist=gamlist,d=3,tol=Tol)$integral)   
  ################ T4
  T4<-((IBR(anr1) - I2BR(1,anr1))*(1 - 2*I2BR(1,anr2)) + (IBR(anr2) -  I2BR(1,anr2))*(1 - 2*I2BR(1,anr1)) -
         adaptIntegrate(I3,lowerLimit=c(0,0),upperLimit=c(1,1),coord=coord,tau=tau,alpha=alpha,coordt=coord[,c(1,3)])$integral - 
         adaptIntegrate(I3,lowerLimit=c(0,0),upperLimit=c(1,1),coord=coord,tau=tau,alpha=alpha,coordt=coord[,c(2,4)])$integral)
  ################# T5
  T5<-((IBR(anr1) - I2BR(1,anr1))*(1 - 2*I2BR(1,anr2)) + (IBR(anr2) -  I2BR(1,anr2))*(1 - 2*I2BR(1,anr1)) -
         adaptIntegrate(I3,lowerLimit=c(0,0),upperLimit=c(1,1),coord=coord,tau=tau,alpha=alpha,coordt=coord[,c(1,4)])$integral - 
         adaptIntegrate(I3,lowerLimit=c(0,0),upperLimit=c(1,1),coord=coord,tau=tau,alpha=alpha,coordt=coord[,c(2,3)])$integral)
  return(T1 - T2 - T3 + T4 + T5)
}

## case 4: s neq u, t = v
intBRcase4<-function(coord,tau,alpha,Tol){ #coord=(s,t,u,v) where s,t,u,v, are 2-dim columnvectors
  loc1<-coord[,2]-coord[,1]
  loc2<-coord[,4]-coord[,3]
  anr1 <- sqrt(2)*((t(loc1) %*% tau %*% loc1)^(alpha/4))
  anr2 <- sqrt(2)*((t(loc2) %*% tau %*% loc2)^(alpha/4))
  ########### T1
  d<-3
  rlist<-gamlist<-vector('list',length=d)
  coordt<-coord[,-4]
  for(i in 1:d){
    rlist[[i]]<-corrfunc(tau,alpha,coordt,d,i)
    gamlist[[i]]<-sapply(c(1:d)[-i], function(j) sqrt(2*gfunc(tau,alpha,coordt[,i]-coordt[,j])))
  }
  T1<-(IBR(anr1)+IBR(anr2)-adaptIntegrate(ellBR4dv2,lowerLimit=c(0,0,0),upperLimit=c(1,1,1),rlist=rlist,gamlist=gamlist,d=3,ind=2,tol=Tol)$integral)
  ############# T2
  T2<-(2*IBR(anr1)*(I2BR(1,anr2) - 0.5) + 2*I2BR(1,anr2) - 2*IBR(anr2) - 
         adaptIntegrate(I1,lowerLimit=c(0,0,0),upperLimit=c(1,1,1),a=anr2,rlist=rlist,gamlist=gamlist,d=3,ind=3,tol=Tol)$integral -
         adaptIntegrate(I1v2,lowerLimit=c(0,0,0),upperLimit=c(1,1,1),a=anr2,tau=tau,alpha=alpha,coordt=coord[,c(1,2)],ind=3,ind2=c(2,3,1),tol=Tol)$integral)
  ############# T3
  T3<-(2*IBR(anr2)*(I2BR(1,anr1) - 0.5) + 2*I2BR(1,anr1) - 2*IBR(anr1) - 
         adaptIntegrate(I1v3,lowerLimit=c(0,0,0),upperLimit=c(1,1,1),a=anr1, rlist=rlist,gamlist=gamlist,d=3,tol=Tol)$integral - 
         adaptIntegrate(I1v2,lowerLimit=c(0,0,0),upperLimit=c(1,1,1),a=anr1,tau=tau,alpha=alpha,coordt=coord[,c(2,3)],ind=1,ind2=c(1,3,2),tol=Tol)$integral)    
  ################ T4
  T4<-((IBR(anr1) - I2BR(1,anr1))*(1 - 2*I2BR(1,anr2)) + (IBR(anr2) -  I2BR(1,anr2))*(1 - 2*I2BR(1,anr1)) -
         adaptIntegrate(I3,lowerLimit=c(0,0),upperLimit=c(1,1),coord=coord,tau=tau,alpha=alpha,coordt=coord[,c(1,3)])$integral - 
         adaptIntegrate(I3,lowerLimit=c(0,0),upperLimit=c(1,1),coord=coord,tau=tau,alpha=alpha,coordt=coord[,c(2,4)])$integral)
  ################# T5
  T5<-((IBR(anr1) - I2BR(1,anr1))*(1 - 2*I2BR(1,anr2)) + (IBR(anr2) -  I2BR(1,anr2))*(1 - 2*I2BR(1,anr1)) -
         adaptIntegrate(I3,lowerLimit=c(0,0),upperLimit=c(1,1),coord=coord,tau=tau,alpha=alpha,coordt=coord[,c(1,4)])$integral - 
         adaptIntegrate(I3,lowerLimit=c(0,0),upperLimit=c(1,1),coord=coord,tau=tau,alpha=alpha,coordt=coord[,c(2,3)])$integral)
  return(T1 - T2 - T3 + T4 + T5)
}

## case 5: s = u, t = v
intBRcase5<-function(coord,tau,alpha,Tol){ #coord=(s,t,u,v) where s,t,u,v, are 2-dim columnvectors
  loc<-coord[,2]-coord[,1]
  a <- sqrt(2)*((t(loc) %*% tau %*% loc)^(alpha/4))
  ########### T1
  coordt<-coord[,1:2]
  T1<-(2*IBR(a)  - adaptIntegrate(ellBR4dv4,lowerLimit=c(0,0),upperLimit=c(1,1),coord=coordt,tau=tau,alpha=alpha,tol=Tol)$integral)
  ############# T2
  T2<-(2*IBR(a)*(I2BR(1,a) - 0.5) + 2*I2BR(1,a) - 2*IBR(a) - 
         adaptIntegrate(I1v2,lowerLimit=c(0,0,0),upperLimit=c(1,1,1),a=a,tau=tau,alpha=alpha,
                        coordt=coordt,ind=3,ind2=c(1,3,2),tol=Tol)$integral - 
         adaptIntegrate(I1v2,lowerLimit=c(0,0,0),upperLimit=c(1,1,1),a=a,tau=tau,alpha=alpha,
                        coordt=coordt,ind=3,ind2=c(2,3,1),tol=Tol)$integral)
  ############# T3
  T3<-(2*IBR(a)*(I2BR(1,a) - 0.5) + 2*I2BR(1,a) - 2*IBR(a) - 
         adaptIntegrate(I1v2,lowerLimit=c(0,0,0),upperLimit=c(1,1,1),a=a,tau=tau,alpha=alpha,
                        coordt=coordt,ind=1,ind2=c(1,2,3),tol=Tol)$integral - 
         adaptIntegrate(I1v2,lowerLimit=c(0,0,0),upperLimit=c(1,1,1),a=a,tau=tau,alpha=alpha,
                        coordt=coordt,ind=1,ind2=c(1,3,2),tol=Tol)$integral)    
  ################ T4
  T4<-((IBR(a) - I2BR(1,a))*(1 - 2*I2BR(1,a)) + (IBR(a) -  I2BR(1,a))*(1 - 2*I2BR(1,a)) -
         adaptIntegrate(I3,lowerLimit=c(0,0),upperLimit=c(1,1),coord=coord,tau=tau,alpha=alpha,coordt=coord[,c(1,3)])$integral - 
         adaptIntegrate(I3,lowerLimit=c(0,0),upperLimit=c(1,1),coord=coord,tau=tau,alpha=alpha,coordt=coord[,c(2,4)])$integral)
  ################# T5
  T5<-((IBR(a) - I2BR(1,a))*(1 - 2*I2BR(1,a)) + (IBR(a) -  I2BR(1,a))*(1 - 2*I2BR(1,a)) -
         adaptIntegrate(I3,lowerLimit=c(0,0),upperLimit=c(1,1),coord=coord,tau=tau,alpha=alpha,coordt=coord[,c(1,4)])$integral - 
         adaptIntegrate(I3,lowerLimit=c(0,0),upperLimit=c(1,1),coord=coord,tau=tau,alpha=alpha,coordt=coord[,c(2,3)])$integral)
  return(T1 - T2 - T3 + T4 + T5)
}


psidotBR<-function(tau,alpha,pairs){ #nrow(pairs) = number of pairs, ncol(pairs)=4
  npairs<-nrow(pairs)
  loc<-pairs[,3:4]-pairs[,1:2]
  psi<-matrix(,nrow=npairs,ncol=4)
  deriv<-vector(length=4)
  for(i in 1:npairs){
    a <- sqrt(2)*((t(loc[i,]) %*% tau %*% loc[i,])^(alpha/4))
    temp<-sqrt(2)*(alpha/4)*((t(loc[i,]) %*% tau %*% loc[i,])^(alpha/4 - 1))
    deriv[1]<-(a*log((t(loc[i,]) %*% tau %*% loc[i,])^(1/4)))
    deriv[2]<-temp*(loc[i,1]^2)
    deriv[3]<-temp*(loc[i,2]^2)
    deriv[4]<-temp*(2*loc[i,1]*loc[i,2])
    psi[i,]<-deriv*((1/2)*dnorm(a/2) - (1/2)*exp(a^2)*dnorm(-(3/2)*a) + (2/3)*a*exp(a^2)*pnorm(-(3/2)*a))
  }
  return(psi)
}

AsymVarBR<-function(pairs,tau,alpha,Tol){ 
  npairs<-nrow(pairs)
  BRMat<-matrix(,nrow=npairs,ncol=npairs)
  for(i in 1:npairs){
    for(j in i:npairs){
      temp<-matrix(c(pairs[i,],pairs[j,]),nrow=2)
      if(all(temp[,1]==temp[,4]) && (! all(temp[,2]==temp[,3]))){
        BRMat[i,j]<-intBRcase2(cbind(temp[,3:4],temp[,1:2]),tau,alpha,Tol=Tol)
      }else if(all(temp[,1]==temp[,3]) && all(temp[,2]==temp[,4])){
        BRMat[i,j]<-intBRcase5(temp,tau,alpha,Tol=Tol)
      }else if(all(temp[,2]==temp[,4]) && (! all(temp[,1]==temp[,3]))){
        BRMat[i,j]<-intBRcase4(temp,tau,alpha,Tol=Tol)
      }else if(all(temp[,1]==temp[,3]) && (! all(temp[,2]==temp[,4]))){
        BRMat[i,j]<-intBRcase3(temp,tau,alpha,Tol=Tol)
      }else if(all(temp[,2]==temp[,3]) && (! all(temp[,1]==temp[,4]))){
        BRMat[i,j]<-intBRcase2(temp,tau,alpha,Tol=Tol)
      }else if((! all(temp[,1]==temp[,3])) && (! all(temp[,2]==temp[,4]))){
        BRMat[i,j]<-intBRcase1(temp,tau,alpha,Tol=Tol)
      }
      BRMat[j,i]<-BRMat[i,j]
    }
  }
  return(BRMat)
}


MestimatorBR <- function(x, locations, pairIndices, k,Tol, startingValue, Omega,iterate,covMat) {
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
  tau<-parsBRtoTau(startingValue[2:4])
  alpha<-startingValue[1]
  theta <- thetaPilot <- optim(par = c(alpha,tau), fn = tominimizeBR, method="BFGS", totlist = totL,
                               pairs = pairs,w = Omega,control=list(maxit=10000,reltol=1e-12))$par
  if (iterate) {
    Omega <- solve(AsymVarBR(pairs, rbind(c(theta[2],theta[4]),c(theta[4],theta[3])), theta[1], Tol))
    theta <- optim(par = thetaPilot,fn = tominimizeBR,method="BFGS",totlist = totL,
                   pairs = pairs,w = Omega,control=list(maxit=10000,reltol=1e-12))$par    
  }
  covMatrix<-NULL
  if(covMat){
    taum<-rbind(c(theta[2],theta[4]),c(theta[4],theta[3]))
    psidot <- psidotBR(taum,theta[1], pairs)
    temp <- solve(t(psidot) %*% Omega %*% psidot) %*% t(psidot)
    GammaMat <- AsymVarBR(pairs, taum, theta[1], Tol)
    covMatrix <- (temp %*% Omega %*% GammaMat %*% Omega %*% t(temp)) / k
  }
  return(list(theta = c(theta[1],parsTauToBR(theta[2:4])),
              theta_pilot = c(thetaPilot[1],parsTauToBR(thetaPilot[2:4])),
              covMatrix = covMatrix,
              Omega = Omega))
}

