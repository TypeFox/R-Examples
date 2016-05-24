## dcl.estimation provide the estimated DCL parameters from two triangles (N,X)
dcl.estimation<-function(Xtriangle,Ntriangle,adj=1,Tables=TRUE,num.dec=4,
                         n.cal=NA,Fj.X=NA,Fj.N=NA) 
{
  Xtriangle<-as.matrix(Xtriangle)
  Ntriangle<-as.matrix(Ntriangle)
  m<-nrow(Xtriangle)
  # m is the dimension of the triangle and d=m-1 is the maximum delay
  d<-m-1
  
  if (nrow(Ntriangle)!=m) stop("The triangles have different dimmension")
  
  ## 1. Estimate the cl parameters
  clm.N<-clm(Ntriangle,n.cal=n.cal,Fj=Fj.N)
  alpha.N<-clm.N$alpha
  beta.N<-clm.N$beta
  Nhat<-as.matrix(clm.N$triangle.hat)
  
  clm.X<-clm(Xtriangle,n.cal=n.cal,Fj=Fj.X) 
  alpha.X<-clm.X$alpha
  beta.X<-clm.X$beta
  
  Xhat<-as.matrix(clm.X$triangle.hat)
  
  ## 2. Delay parameters estimation (by solving the linear-equations system (7) 
  #     in the paper and removing possible NA's
  nonas<-(!is.na(beta.X) ) & (!is.na(beta.N)) 
  beta.X<-beta.X[nonas]
  beta.N<-beta.N[nonas]
  mm<-length(beta.X)
  ## create the matrix of beta's in the system'
  mat.beta.N<-matrix(0,mm,mm)
  mat.beta.N[,1]<-beta.N
  gen.mat<-function(l) mat.beta.N[,l]<-c(rep(0,length((mm-l+2):mm)),beta.N[-((mm-l+2):mm)])
  mat.beta.N[,2:mm]<-apply(t(2:mm),2,gen.mat)
  
  ## solving the system with the function 'solve' doesn't always work
  ## we calculate the general inverse
  g.inv<-function(X, tol = sqrt(.Machine$double.eps))
  {
    ## Generalized Inverse of a Matrix
    dnx <- dimnames(X)
    if(is.null(dnx)) dnx <- vector("list", 2)
    s <- svd(X)
    nz <- s$d > tol * s$d[1]
    structure(
      if(any(nz)) s$v[, nz] %*% (t(s$u[, nz])/s$d[nz]) else X,
      dimnames = dnx[2:1])
  }
  betainv<-g.inv(mat.beta.N)
  p.delay<-betainv%*%beta.X
  
  ## warning if the DCL model is not working
  if (max(abs(p.delay))>1) {
    message("Warning: It is not possible to extract a suitable delay function. Check if the data follow the DCL model")
  }
  ## now we adjust the solution because 'p.delay' could be not a probability vector
  
  # this is the method used for the illustration in the DCL paper
  if (adj==1) {
    if (any(p.delay<0))
    {
      ind.pos<-min(which(p.delay<0),na.rm=T)
    } else ind.pos<-d+1 
    ## then hold all the elements but normalize to add up 1 with the last
    p.delay.DCL<-p.delay[1:(ind.pos-1)]
    Sp.delay<-cumsum(p.delay.DCL)
    ind.pos<-(Sp.delay<1) 
    ## we reduce the maximum number of periods delay to d= length(ind.pos)+1
    p.delay.DCL<-p.delay.DCL[ind.pos]
    p.delay.DCL<-c(p.delay.DCL,1-sum(p.delay.DCL))  
    pj<-p.delay.DCL
    d.dcl<-length(pj)
    ## get the maximum delay d=m-1 
    pj<-c(pj,rep(0,m-d.dcl))
  }
  
  # other way is proportional to the size of the pj's
  # it works better with DCL in this case
  if (adj==2) {
    pj<-p.delay
    pj[p.delay<0]<-0
    rest<-(1-sum(pj))
    v.rest<-rest*pj/sum(pj);pj<-pj+v.rest
  }
  
  #### 3. Inflation and severity parameters estimation
  
  # The following was added to work with triangles having first row with zeros
  ind.pos<-which(alpha.X!=0 & alpha.N!=0)
  set.one<-min(ind.pos)
  ## Let inf.hat[set.one]=1 then
  mu<- alpha.X[set.one]/alpha.N[set.one]
  ## the estimated inflation parameter
  inflat<-alpha.X/(mu*alpha.N)
  if (set.one>1) inflat[1:(set.one-1)]<-0
  # if we model the delay parameters as a probability vector (pj)
  # to ensure that the beta's sum 1 (identification) we need to adjust the mean by constant
  new.beta.X<-mat.beta.N%*%pj
  kappa<-sum(new.beta.X)
  mu.adj<-mu/kappa
  
  ## We will have problems if the alpha.X is zero in any row so we include a filter 
  for (i in 2:m) if (is.na(inflat[i])) inflat[i]<-inflat[i-1]
  #  with the above sentence we avoid errors when some alpha.N=0 or some alpha.X=0
  
  ## Estimate the variance using the overdispersion property
  #  same as VNJ (2010) once the data has been normalized removing the inflation
  ## We will have problems if the alpha.X is zero in any row so we include a filter 
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
  
  ## the estimated means (different at each accident year)
  Ey<-mu.adj *inflat #inf.corr
  
  if (Tables==TRUE)
  {  
    ## print some tables
    table1<-data.frame(delay.par=p.delay, 
                       delay.prob=pj,inflation=inflat,
                       severity.mean=Ey,severity.var=Vy )
    table2<-data.frame(mean.factor=mu,mean.factor.adj=mu.adj,
                       variance.factor=sigma2)
    ## Print the results
    print(round(table1,num.dec))
    print(round(table2,num.dec))
  }
  
#   if (Plots==TRUE)
#   {
#     g.dev
#     par(mfrow=c(2,2))
#     par(mar=c(3,1.5,1.5,1.5),oma=c(2,0.5,0.5,0.2),mgp=c(1.5,0.5,0))
#     # c(bottom, left, top, right)
#     
#     par(cex.main=1)
#     plot(1:m,alpha.X,type='b' ,pch=19,col=1,cex.axis=0.8,xlab='underwriting period',
#          ylab='', main='CL underwriting parameters',xaxt='n')
#     axis(1,1:m,as.character(1:m))
#     
#     grid(ny=5,nx=NA)
#     plot(1:m,beta.X,type='b' ,pch=19,lty=1,col=1,cex.axis=0.8,xaxt='n',xlab='development period',
#          ylab='', main='CL development parameters')
#     axis(1,1:m,as.character(0:(m-1)))
#     grid(ny=5,nx=NA)
#     plot(1:m,inflat,col=4,pch=19,ylab='',xlab='underwriting period',
#          main='Severity inflation',cex.axis=0.8,xaxt='n')
#     axis(1,1:m,as.character(1:m))
#     
#     grid(ny=5,nx=NA)
#     lines(inflat,col=4,lty=2,lwd=1)
#     plot(1:m,p.delay,type='l',lwd=2,lty=1,col=4,xaxt='n',xlab='settlement delay',
#          ylab='',main='Delay parameters',cex.axis=0.8)
#     axis(1,1:m,as.character(0:(m-1)))
#     grid(ny=5,nx=NA)
#     points(1:m,p.delay,pch=18,col=4,cex=0.6)
#     lines(1:m,pj,type='l',lwd=2,lty=2,col=3)
#     points(1:m,pj,pch=15,col=3,cex=0.6)
#     grid(ny=5,nx=NA)
#     legend('topright',legend=c('general','adjusted'),col=c('blue','green'),
#            bty='n',pch=c(18,15),cex=0.9)#,lty=c(1,2))
#     par(mfrow=c(1,1))
#     
#   }  
  
  # It returns a list with:
  return(list(pi.delay=p.delay,mu=mu,inflat=inflat,pj=pj,mu.adj=mu.adj,
              sigma2=sigma2,phi=phi,Ey=Ey,Vy=Vy, adj=adj,
              alpha.N=alpha.N,beta.N=beta.N,Nhat=Nhat,
              alpha.X=alpha.X,beta.X=beta.X,Xhat=Xhat))
}