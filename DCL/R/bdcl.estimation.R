bdcl.estimation<-function(Xtriangle,Ntriangle,Itriangle,adj=1,
                          Tables=TRUE,num.dec=4,n.cal=NA,Fj.X=NA,Fj.N=NA,Fj.I=NA) 
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
  Xhat<-as.matrix(par.dcl$Xhat)
  alpha.X<-par.dcl$alpha.X
  beta.X<-par.dcl$beta.X
  inflat.DCL<-par.dcl$inflat
  alpha.N<-par.dcl$alpha.N
  beta.N<-par.dcl$beta.N
  Nhat<-as.matrix(par.dcl$Nhat)
  
  ## Get the inflation from incurred data
  clm.I<-clm(Itriangle,n.cal=n.cal,Fj=Fj.I)  
  alpha.I<-clm.I$alpha
  beta.I<-clm.I$beta
  
  #mu.I<-alpha.I[1]/alpha.N[1] #it is the same as mu
  
  # filter to work with triangles which have first row with all zeros
  ind.pos<-which(inflat.DCL!=0)
  set.one<-min(ind.pos)
  ## Let inf.hat[set.one]=1 then
  ## the estimated inflation parameter
  inflat<-alpha.I/(mu*alpha.N)
  if (set.one>1) inflat[1:(set.one-1)]<-0
  
  ## We will have problems if the alpha.X is zero in any row so we include a filter 
  for (i in (set.one:m)) if (is.na(inflat[i])) inflat[i]<-inflat[i-1]
  #  with the above sentence we avoid errors when some alpha.N=0
  
  ## Estimate the variance using the overdispersion property
  #  same as VNJ (2010) once the data has been normalized removing the inflation
  ## We will have problems if the alpha.X is zero in any row so we include a filter 
  sigma2DCL<-par.dcl$sigma2
  
  if (!is.na(sigma2DCL))
  {  
    inf.corr<-inflat
    inf.corr[inf.corr==0]<-1 
    Ztriangle<-Xtriangle
    for (i in 1:m) Ztriangle[i,]<-Xtriangle[i,]/inf.corr[i]
    
    ## the conditional expectation (ignoring inflation i.e. inflat=rep(1,m) ) 
    Ec.setI<-function(Ntriangle,pj,mu,inflat,m,d) 
    {
      v.expect<-apply(expand.grid(1:m,1:m),MARGIN=1,
                      FUN= function(v)
                      { j<-v[1];i<-v[2];
                        if (j<(m-i+2)) v.e<-sum(Ntriangle[i,((max(j-d,1)):j)] 
                                                * pj[((j+1)-((max(j-d,1)):j))] * mu*inflat[i]) else v.e<-0
                        return(v.e) 
                      })
      expect<-matrix(v.expect,m,m,byrow=T)
      # It returns the triangle of dimension m with the conditional expectation 
      return(expect)
    }
    expect<-Ec.setI(Ntriangle,pj,mu.adj,rep(1,m),m,d) 
    ## the number degree-of-freedom (maybe some terms are zero)
    df_phi<-sum(expect!=0)-(d+1)
    phi<-sum((Ztriangle[expect!=0]-expect[expect!=0])^2/expect[expect!=0],na.rm=T)/df_phi
    sigma2<- mu.adj*(phi- mu.adj)
    if (is.na(sigma2)){
      message('There is not enough information to estimate the variance of the 
              individual size claims. See the documentation for suggestions.')
      Vy<-NA;sigma2<-NA
    } else {
      if (sigma2<=0){
        message('There is not enough information to estimate the variance of the 
              individual size claims. See the documentation for suggestions.')
        Vy<-NA;sigma2<-NA
      } 
      Vy<-inf.corr^2 *sigma2
    }
  } else {Vy<-NA;sigma2<-NA}
  
  ## the estimated means (different at each accident year)
  Ey<-mu.adj*inflat #inf.corr
  
  if (Tables==TRUE)
  {  
    ## print some tables
    table1<-data.frame(delay.par=pi.delay, 
                       delay.prob=pj,inflat.DCL=inflat.DCL,inflat.BDCL=inflat,
                       severity.mean=Ey,severity.var=Vy )
    table2<-data.frame(mean.factor=mu,mean.factor.adj=mu.adj,
                       variance.factor=sigma2)
    ## Print the results
    print(round(table1,num.dec))
    print(round(table2,num.dec))
  }

  # It returns a list with:
  return(list(pi.delay=pi.delay,mu=mu,inflat=inflat,inflat.DCL=inflat.DCL,pj=pj,mu.adj=mu.adj,
              sigma2=sigma2,phi=phi,Ey=Ey,Vy=Vy, adj=adj,
              alpha.N=alpha.N,beta.N=beta.N,Nhat=Nhat,
              alpha.X=alpha.X,beta.X=beta.X,Xhat=Xhat,
              alpha.I=alpha.I,beta.I=beta.I))
}