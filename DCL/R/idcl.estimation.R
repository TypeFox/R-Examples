idcl.estimation<-function(Xtriangle,Ntriangle,Itriangle,adj=1,Tables=TRUE,
                          num.dec=4,n.cal=NA,Fj.X=NA,Fj.N=NA,Fj.I=NA)  
{
  Ntriangle<-as.matrix(Ntriangle)
  Xtriangle<-as.matrix(Xtriangle)
  Itriangle<-as.matrix(Itriangle)
  m<-nrow(Ntriangle)
  ## DCL parameters from triangles X (paid) and N -original DCL 
  par.dcl<-dcl.estimation(Xtriangle,Ntriangle,adj,Tables=FALSE,n.cal=n.cal,
                          Fj.X=Fj.X,Fj.N=Fj.N)
  pj<-par.dcl$pj
  pi.delay<-par.dcl$pi.delay
  d<-length(pj)-1
  mu<-par.dcl$mu
  mu.adj<-par.dcl$mu.adj
  inflat.DCL<-par.dcl$inflat
  alpha.N<-par.dcl$alpha.N
  beta.N<-par.dcl$beta.N
  Nhat<-as.matrix(par.dcl$Nhat)
  alpha.X<-par.dcl$alpha.X
  beta.X<-par.dcl$beta.X
  Xhat<-par.dcl$Xhat
  low.Xhat<-as.matrix(Xhat)
  low.Xhat[row(Xhat)+col(Xhat)<=(m+1)]<-0
  CL.X.i<-rowSums(low.Xhat,na.rm=T)
  
  ## Get the inflation from incurred data to reproduce exaclty incurred reserve
  clm.I<-clm(Itriangle,n.cal=n.cal,Fj=Fj.I) 
  alpha.I<-clm.I$alpha
  beta.I<-clm.I$beta
  # The total paid for each accident year in the past
  Ri.X<-rowSums(Xtriangle,na.rm=T)
  # Incurred outstanding numbers
  CL.I.i<-alpha.I-Ri.X
    
  inflat.factor<-CL.I.i/CL.X.i
  inflat.factor[CL.X.i==0]<-1
  inflat<-inflat.factor*inflat.DCL
  inflat[abs(inflat)==Inf]<-NA
  
  ## Do not estimate again the variance, use the DCL
  inf.corr<-inflat
  inf.corr[inf.corr==0]<-1 
  sigma2<-par.dcl$sigma2
  phi<-par.dcl$phi
  
  Vy<-inf.corr^2 *sigma2
  ## the estimated means (different at each accident year)
  Ey<-mu.adj*inflat #inf.corr
  
  if (Tables==TRUE)
  {  
    ## print some tables
    table1<-data.frame(delay.par=pi.delay, 
                       delay.prob=pj,inflat.DCL=inflat.DCL,inflat.IDCL=inflat,
                       severity.mean=Ey,severity.var=Vy )
    table2<-data.frame(mean.factor=mu,mean.factor.adj=mu.adj,
                       variance.factor=sigma2)
    ## Print the results
    print(table1)
    print(table2)
  }
    
  # It returns a list with:
  return(list(pi.delay=pi.delay,mu=mu,inflat=inflat,inflat.DCL=inflat.DCL,pj=pj,mu.adj=mu.adj,
              sigma2=sigma2,phi=phi,Ey=Ey,Vy=Vy, adj=adj,
              alpha.N=alpha.N,beta.N=beta.N,Nhat=Nhat,
              alpha.X=alpha.X,beta.X=beta.X,Xhat=Xhat,
              alpha.I=alpha.I,beta.I=beta.I,CL.I.i=CL.I.i))
}
