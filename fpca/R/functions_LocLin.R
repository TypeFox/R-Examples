###functions for local linear fitting
###5-22-07

##########################get inital estimate of B, lambda, sig by pooled local linear method
###format data.list in terms of Obs, T and N
Format.data<-function(data.list,n,Nmax){
##para: data.list: result of "GenObs"
##n: number of subjects; Nmax: max number of measurements per subject
Obs<-matrix(0,n,Nmax)
T<-matrix(0,n,Nmax)
N<-numeric(Nmax)
  for(i in 1:n){
   cur<-matrix(data.list[[i]][[1]],ncol=2)
   N[i]<-nrow(cur)
   Obs[i,1:N[i]]<-cur[,1]
   T[i,1:N[i]]<-cur[,2]
  }
  return(list(Obs,T,N))
}

###
TranMtoV<-function(Obs,L,N){
##transform data into a vector
##name:TranMtoV 
##para: Obs--observed values,L--locations of measurements, N--# of measurements of each curve
##return: data--3col matrix;1--ID,2--data,3--locations of meaturements
 n<-sum(N)
 result<-matrix(0,n,3)
 count<-0
 for(i in 1:nrow(Obs)){
  for(j in 1:N[i]){
   count<-count+1
   result[count,1]<-i
   result[count,2]<-Obs[i,j]
   result[count,3]<-L[i,j]
   }
  }
return(result)
}



##
###
LocLinEst.new<-function(Obs,T,N,indext,hmu,hcov,hsig,epsi=1e-8){
##local linear estimation of mu(t), G(s,t) for given bandwidth, and evaluated at fine grids
##name:LocLinEst
##para:Obs--observed values, T--measurement time points, N--number of measurements;
##     indext--1 dim evaluation grid;
##     hmu--bandwidth for mu,hcov--bandwidth for G(s,t)(2*1),hd--bandwidth for diagonal,
##     hsig--bandwidth for G(t,t)+sigma^2;
##result:ssmt$estimate--estimation of mu; covmatrix--estimation of G(s,t)I(s!=t);
##       diagcovt--estimation of G(t,t) with quadratic kernel by CovDiag; 
##       ssmdt$estimate--estimation of G(t,t)+\sigma^2 by LocLinErr.diag             
##       sighat2--estimation of sigma^2*(L2-L1); 

  data<-TranMtoV(Obs,T,N)  #transform data to vectors  
  trdata<-TransData(Obs,T,N)
  y<-data[,2]              #obs
  t<-data[,3]              #measurements 
#(0) two dim evaluation grid
  t.points<-indext    #2-dim
  s.points<-indext
  M<-length(t.points)
  x.points<-matrix(0,2,M*(M+1)/2)
  count<-0
  for (i in 1:M){
   for (j in i:M){
    count<-count+1
    x.points[1,count]<-t.points[i]
    x.points[2,count]<-s.points[j]
   }
  }

#(i) estimate mean by local linear smoother:
  ssm<-sm.regression(t,y,h=hmu,poly.index=1,eval.points=t,eval.grid=FALSE)       #evaluate at t
  ssmt<-sm.regression(t,y,h=hmu,poly.index=1,eval.points=indext,eval.grid=FALSE) 
  #evaluate at equally spaced grid
  ssmtr<-sm.regression(t,y,h=hmu,poly.index=1,eval.points=as.vector(trdata[[3]]),eval.grid=FALSE) 

#(iii)estimate covariance: G(s,t)I(s!=t)
  fitmu<-ssm$estimate   
  conver<-HatGOff(fitmu,data,N)
  Cova<-conver[[2]]
  CovT<-conver[[3]]  
  ssmcovt<-sm.regression(t(CovT),Cova,h=hcov,poly.index=1,eval.points=t(x.points),eval.grid=FALSE) 

#(iv)covert to matrix
  covmatrix<-matrix(0,M,M)  #G(s,t)I(s!=t) (linear diagonal) 
  count<-0
  for (i in 1:M){
    for (j in i:M){
     count<-count+1
     covmatrix[i,j]<-ssmcovt$estimate[count]
     covmatrix[j,i]<-covmatrix[i,j]           
    }    
  }
  
#(v) estimate G(t,t)+\sigma^2 and \sigma^2   
 ssmdt<-sm.regression(t,(y-fitmu)^2,h=hsig,poly.index=1,eval.points=indext,eval.grid=FALSE) 
 diagsig<-ssmdt$estimate

 sig2hat<-sum(diagsig[round(M/4):round(3*M/4)]-diag(covmatrix)[round(M/4):round(3*M/4)])*(indext[2]-indext[1])*2

#(vi) return 
  result<-list(ssmt$estimate,covmatrix,diagsig,sig2hat)
  return(result)
  }



###
HatGOff<-function(fitmu,data,N){
##standardize data by subtracting mean and then return \hat{G(s,t)} at off diagonal pair (s,t)
##result used by LocLin.G; called by CV.G
##name:HatGOff
##para:fitmu--fitted mean curve at observed points; 
##data<-TranMtoV(Obs,T,N) format;N--number of observations for each subject
##return:Cova--a vector of \hat{G(tij,til)}; and CovT: time pairs of observation (tij,til),j<l; 
##       index--subject ID for the pair (tij,til), j<l 
  data.mean<-cbind(data,fitmu)
  Cova<-numeric(sum(N*(N-1)/2))     
  CovT<-matrix(0,nrow=2,ncol=sum(N*(N-1)/2))                                    
  index<-numeric(sum(N*(N-1)/2)) 
   count<-0  
    for (i in 1:length(N)){
      if(N[i]>1){
      med<-data.mean[data[,1]==i,]
      yi<-med[,2]
      ti<-med[,3] 
      mi<-med[,4]
       for (j in 1: (N[i]-1)){
        for (l in (j+1): N[i]){
         count<-count+1
          Cova[count]<-(yi[j]-mi[j])*(yi[l]-mi[l])
          CovT[1,count]<-ti[j] 
          CovT[2,count]<-ti[l] 
          index[count]<-i
         }
       }
     }
    }  
  return(list(index,Cova,CovT))
}

###
TransData<-function(Obs,L,N){ 
##transform data to 2 column matrix (Y_{ij},Y_{il},j!=l)and(L_{ij},L_{il},j<l);
##then apply diagonal transformation on Ls;
##result used by CovDiag; called by CV.diag
##Name:TransData
##para: Obs--observed values,L--locations of measurement, N--# of measurements of each curve,
##return: TObs--observation(2 columns), TL--tranformed times;
##        TLO--original times, TId--subject id for time each pair
  tran<-sqrt(2)/2*matrix(c(1,1,-1,1),2,2,byrow=T)
  num<-sum(N*(N-1)/2) 
  TObs<-matrix(0,num,2)
  TL<-matrix(0,num,2)
  TLO<-matrix(0,num,2)
  TId<-numeric(num)
  count<-0
 for (i in 1:nrow(Obs)){
   for(j in 1:N[i]){
    for(l in j:N[i]){
      if(l !=j){
          count<-count+1   
          TLO[count,]<-c(L[i,j],L[i,l])     
          TL[count,]<-tran%*%matrix(TLO[count,],2,1)                      
          TObs[count,]<-c(Obs[i,j],Obs[i,l])
          TId[count]<-i     
         }     
       }
     }  
   }
return(list(TObs,TL, TLO,TId))
}

##local linear as initial value
LocLin.Ini<-function(data.list,n,nmax,grid.l,grids,r){

###formatting
 data.f<-Format.data(data.list,n,Nmax=nmax)
 Obs<-data.f[[1]]
 T<-data.f[[2]]
 N<-data.f[[3]]
 data<-TranMtoV(Obs,T,N)
 y<-data[,2]
 t<-data[,3]

###bandwidth selection
#(i) mu
 hmu.cv<-h.select(t,y,method="cv")
 fitmu<-sm.regression(t,y,h=hmu.cv,poly.index=1,eval.points=t)$estimate 
# plot(t,fitmu,ylim=range(y),xlab="t",ylab="fitted mu(t)",main="estimated mean curve")
# points(t,y,type='p')

#G(s,t)I(s!=t)
 conver<-HatGOff(fitmu,data,N)  ##convert to off diagonal pair forma
 Cova<-conver[[2]]
 CovT<-conver[[3]]
 hcov.cv<-h.select(t(CovT),Cova,method="cv") 

#G(t,t)+\sigma^2
hsig.cv<-h.select(t,(y-fitmu)^2,method="cv")

###local linear fitting
lin.result<-LocLinEst.new(Obs,T,N,grid.l,hmu.cv,hcov=hcov.cv,hsig.cv)

###result 
 muhat<-lin.result[[1]]
 covmatrix<-lin.result[[2]]
 sig2hat<-lin.result[[4]]/(1-0)

###project on a finer grid
 timeindex<-floor(grids*length(grid.l))+1
 timeindex[length(grids)]<-length(grid.l)
 covmatrix<-covmatrix[timeindex,timeindex]

##get eigenfunctions and eigenvalues 
 eigen.c<-EigenC(covmatrix,grids)
 eigenf.c<-eigen.c[[1]][,1:r]
 eigenv.c<-eigen.c[[2]][1:r]

##selected bandwidth
 band<-list("h.mu"=hmu.cv,"hcov"=hcov.cv,"hsig"=hsig.cv)

##
 return(list(sig2hat,eigenf.c,eigenv.c,band))
}


###
#local linear to estimate the mean; then substract it from data,
#then return data.list (mean subtracted)
LocLin.mean<-function(data.list,n,nmax,grids=seq(0,1,0.002)){

###formatting
 data.f<-Format.data(data.list,n,Nmax=nmax)
 Obs<-data.f[[1]]
 T<-data.f[[2]]
 N<-data.f[[3]]
 data<-TranMtoV(Obs,T,N)
 y<-data[,2]
 t<-data[,3]

###bandwidth selection
#(i) mu
 hmu.cv<-h.select(t,y,method="cv")
 fitmu<-sm.regression(t,y,h=hmu.cv,poly.index=1,eval.points=t)$estimate 
 fitmu.s<-sm.regression(t,y,h=hmu.cv,poly.index=1,eval.points=grids)$estimate 

##subtract mean 
 data.s<-data
 data.s[,2]<-data[,2]-fitmu  

##write back
 data.result<-data.list 
 for(i in 1:n){
  data.c<-matrix(data.s[data.s[,1]==i,],ncol=3)
  data.result[[i]][[1]]<-data.c[,2:3]
 }

 return(list(data.result,fitmu.s))
}
