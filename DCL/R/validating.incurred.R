## The following function performs a backtest in the considered aggregation period
validating.incurred<-function(ncut=0,Xtriangle,Ntriangle,Itriangle,Model=0,
                              Plot.box=TRUE,Tables=TRUE,num.dec=4, 
                              n.cal=NA,Fj.X=NA,Fj.N=NA,Fj.I=NA)
{ 
  Ntriangle<-as.matrix(Ntriangle)
  Xtriangle<-as.matrix(Xtriangle)
  Itriangle<-as.matrix(Itriangle)
  m<-nrow(Ntriangle)
  ncut<-as.integer(ncut)
  if (ncut>m | ncut<0) stop("Value ncut not allowed")
  if (ncut==0)
  {  ## Estimate the three inflations from the full triangle
    clm.N<-clm(Ntriangle,n.cal=n.cal,Fj=Fj.N)
    alpha.N<-clm.N$alpha
    beta.N<-clm.N$beta
    
    clm.X<-clm(Xtriangle,n.cal=n.cal,Fj=Fj.X) 
    alpha.X<-clm.X$alpha
    beta.X<-clm.X$beta
    Xhat<-as.matrix(clm.X$triangle.hat)
    
    
    clm.I<-clm(Itriangle,n.cal=n.cal,Fj=Fj.I)
    Ihat<-clm.I$triangle.hat
    alpha.I<-clm.I$alpha
    
    ## DCL inflation
    # The following was added to work with triangles having first row with zeros
    ind.pos<-which(alpha.X!=0 & alpha.N!=0)
    set.one<-min(ind.pos)
    ## Let inf.hat[set.one]=1 then
    mu<- alpha.X[set.one]/alpha.N[set.one]
    ## the estimated inflation parameter
    inflat.DCL<-alpha.X/(mu*alpha.N)
    #if (set.one>1) inflat.DCL[1:(set.one-1)]<-0
    
    ## BDCL inflation
    inflat.BDCL<-alpha.I/(mu*alpha.N)
    #if (set.one>1) inflat.BDCL[1:(set.one-1)]<-0
    
    ## IDCL inflation
    # The total paid for each accident year in the past
    Ri.X<-rowSums(Xtriangle,na.rm=T)
    # Incurred outstanding numbers
    CL.I.i<-alpha.I-Ri.X
    low.Xhat<-Xhat
    low.Xhat[row(Xhat)+col(Xhat)<=(m+1)]<-0
    CL.X.i<-rowSums(low.Xhat,na.rm=T)
    inflat.factor<-CL.I.i/CL.X.i
    inflat.factor[abs(CL.X.i)<1e-6]<-1   ## this is the way in which RSA does CL.incurred
    inflat.IDCL<-inflat.factor*inflat.DCL
        
    ## Plot with the three inflations
    #if(dev.cur() == 1) dev.new()
    par(mfrow=c(2,1))
    par(mar=c(3,1.5,2.0,1.5),oma=c(2,0.5,0.5,0.2),mgp=c(1.5,0.5,0))
    # # c(bottom, left, top, right)
    par(cex.axis=0.9,cex.main=1)
   
    yy<-range(c(inflat.DCL,inflat.IDCL,inflat.BDCL),na.rm=T)
    plot(1:m,inflat.DCL,col=2,type='b',pch=19,ylab='',xlab='underwriting period',
         main='Severity inflation',xaxt='n',ylim=yy,lwd=2)
    axis(1:m,as.character(1:m))
    grid(ny=5,nx=NA)
    points(inflat.BDCL,col=4,pch=15)
    lines(inflat.BDCL,col=4,lty=1,lwd=2)
    lines(inflat.IDCL,col=3,lty=3,lwd=2)
    points(inflat.IDCL,col=3,pch=18)
    grid(ny=5,nx=NA)
    legend('topleft',c('DCL','BDCL','IDCL'),
           col=c(2,4,3),pch=c(19,15,18),cex=0.7,bty='n')
    
    plot(1:m,inflat.factor,type='b',lwd=1,lty=1,col=1,xlab='underwriting period',
         ylab='',main='Correction factor (IDCL/DCL)',xaxt='n')
    axis(1:m,as.character(1:m))
    grid(ny=5,nx=NA)
    points(1:m,inflat.factor,pch=18,col=1)
    #legend('bottomright','IDCL/DCL)',col=1,
    #       cex=0.7,bty='n',pch=18)
    par(mfrow=c(1,1))
    return(list(inflat.DCL=inflat.DCL,
                inflat.BDCL=inflat.BDCL,inflat.IDCL=inflat.IDCL))
  } else {
  
    ## now we do the backtest, start by creating the subtriangle for the estimation
    cut.diagonals<-function(my.triangle,m)
    {
      # m is the desired dimension of the cutted triangle  
      triangle<-matrix(NA,m,m)
      for (i in 1:m)
      {
        for (j in 1:m)
        {  if (i+j-1<=m) triangle[i,j]<-my.triangle[i,j] else triangle[i,j]<-NA
        }
      }
      triangle
    }
    #get.calendar produces vector of diagonal sums
    get.calendar<-function(triangle,m,ncut)
    {
      calendarvec<-vector(length=ncut)
      for (i in 1:ncut)
      {
        calendarvec[i]<-sum(triangle[row(triangle)+col(triangle)==m+i+1])
      }
      return(calendarvec)
    }
    
    #gives vector of the diagonals from m+1 to ncut
    get.diagonalcuts<-function(triangle,m,ncut)
    {
      diagvec<-triangle[row(triangle)+col(triangle)==m+2]
      if (ncut>1)
        for (i in 2:ncut)
          diagvec<-c(diagvec,triangle[row(triangle)+col(triangle)==m+i+1])
      return(diagvec)
    }
    
    ## The dimensions
    mm<-m  ## the original dimension
    m<-mm-ncut  ## the dimension of the subtriangle from which the estimation is performed
    d<-m-1
    
    ## a) First keep the observed actual values from Xtriangle
    # work with the cutted data but keep the true cells where we will compare
    # the projections by different methods
    chop.obs.X<- Xtriangle
    chop.obs.X[row(Xtriangle)+col(Xtriangle)<=(m+1)]<-NA  
    chop.obs.X[row(Xtriangle)+col(Xtriangle)>(mm+1)]<-NA
    # we cannot project from the row mth so we remove these rows in the comparisons 
    chop.obs.X<-chop.obs.X[1:m,]  # dim =(m,mm)
    
      
    ## b) Now cut the last ncut calendar periods
    Ntriangle<-cut.diagonals(Ntriangle,m)
    Xtriangle<-cut.diagonals(Xtriangle,m)
    Itriangle<-cut.diagonals(Itriangle,m)
    
    
    # 1) project now using DCL with N and X
    
    par.N<-clm(Ntriangle,Fj=Fj.N)  
    #it returns: list(triangle.hat,alpha,beta,clm)
    Nhat<-as.matrix(par.N$triangle.hat)
    alpha.N<-par.N$alpha
    beta.N<-par.N$beta
    
    par.estim.DCL<-dcl.estimation(Xtriangle,Ntriangle,adj=1,Tables=FALSE,
                                  n.cal=n.cal,Fj.X=Fj.X,Fj.N=Fj.N)
    inflat.DCL<-par.estim.DCL$inflat
    mu<-par.estim.DCL$mu
    Xhat<-as.matrix(par.estim.DCL$Xhat)
    alpha.X<-as.vector(par.estim.DCL$alpha.X)
    Ri.X<-rowSums(Xtriangle,na.rm=T)
    CL.X.i<-alpha.X-Ri.X
    
    
    preds<-dcl.predict(par.estim.DCL,Ntriangle,Model=Model,Tail=TRUE,Tables=FALSE)
        
    Xtotal.DCL<-as.matrix(preds$Xtotal)
    Xtotal.DCL.mod<-Xtotal.DCL[,1:m]
    CL.X.i<-rowSums(Xtotal.DCL.mod,na.rm=T)
    
    #Xtotal.DCL[is.na(chop.obs.X)]<-NA
    ## no remove the tail
    ## compare with the actual numbers in chop.obs.X
    #chop.obs.X  #dim =(m,2*m+1)=(14,32)
    chop.obs.X<-cbind(chop.obs.X, matrix(NA,m,ncol(Xtotal.DCL)-ncol(chop.obs.X)))
    chop.obs.X.diag<-get.calendar(chop.obs.X,m,ncut)
    # now dim(chop.obs.X)= (m,m+d)
    
    XdifDCL<-(Xtotal.DCL-chop.obs.X)
    dif.sq<-(XdifDCL)^2
    point.pe.DCL<-sqrt((sum(dif.sq,na.rm=T))/sum(chop.obs.X^2,na.rm=T))
    
    #total.err.DCL<-sum(Xtotal.DCL, na.rm=T)-sum(chop.obs.X, na.rm=T)
    total.err.DCL<-(sum(XdifDCL, na.rm=T))
    total.err.DCL<-abs(total.err.DCL)
    total.pe.DCL<-total.err.DCL/sum(chop.obs.X,na.rm=T)
    
    calendar.dif<-get.calendar(XdifDCL,m,ncut)
    calendar.pe.DCL<-sqrt((sum(calendar.dif^2,na.rm=T))/sum(chop.obs.X.diag^2,na.rm=T))
    
    ## 2) BDCL 
    ## the BDCL inflation is that we did using DCL over N and I
    clm.I<-clm(Itriangle,n.cal=n.cal,Fj=Fj.I) 
    alpha.I<-clm.I$alpha
    # filter to work with triangles which have first row with all zeros
    ind.pos<-which(inflat.DCL!=0)
    set.one<-min(ind.pos)
    ## Let inf.hat[set.one]=1 then
    ## the estimated inflation parameter
    inflat.BDCL<-alpha.I/(mu*alpha.N)
    if (set.one>1) inflat[1:(set.one-1)]<-0
    ## We will have problems if the alpha.X is zero in any row so we include a filter 
    for (i in (set.one:m)) if (is.na(inflat.BDCL[i])) inflat.BDCL[i]<-inflat.BDCL[i-1]
    #  with the above sentence we avoid errors when some alpha.N=0
    
    par.estim.BDCL<-par.estim.DCL
    par.estim.BDCL$inflat<-inflat.BDCL
    preds<-dcl.predict(par.estim.BDCL,Ntriangle,Model=Model,Tail=TRUE,Tables=FALSE)
    
    Xtotal.BDCL<-as.matrix(preds$Xtotal)
    
    #Xtotal.BDCL[is.na(chop.obs.X)]<-NA
    
    XdifBDCL<-(Xtotal.BDCL-chop.obs.X)
    dif.sq<-(XdifBDCL)^2
    point.pe.BDCL<-sqrt((sum(dif.sq,na.rm=T))/sum(chop.obs.X^2,na.rm=T))
    
    #total.err.BDCL<-sum(Xtotal.BDCL, na.rm=T)-sum(chop.obs.X, na.rm=T)
    
    total.err.BDCL<-(sum(XdifBDCL, na.rm=T))
    total.err.BDCL<-abs(total.err.BDCL)
    total.pe.BDCL<-total.err.BDCL/sum(chop.obs.X,na.rm=T)
    
    calendar.dif<-get.calendar(XdifBDCL,m,ncut)
    calendar.pe.BDCL<-sqrt((sum(calendar.dif^2,na.rm=T))/sum(chop.obs.X.diag^2,na.rm=T))
    
    
    # 3) IDCL
    alpha.I<-clm.I$alpha
    # The total paid for each accident year in the past is Ri.X
    # Incurred outstanding numbers
    CL.I.i<-alpha.I-Ri.X
    
    inflat.factor<-CL.I.i/CL.X.i
    inflat.factor[abs(CL.X.i)<1e-6]<-1
    ## The final inflation with IDCL
    inflat.IDCL<-inflat.factor*inflat.DCL
    par.estim.IDCL<-par.estim.DCL
    par.estim.IDCL$inflat<-inflat.IDCL
    preds<-dcl.predict(par.estim.IDCL,Ntriangle,Model=Model,Tail=TRUE,Tables=FALSE)
    
    Xtotal.IDCL<-preds$Xtotal
    #Xtotal.IDCL[is.na(chop.obs.X)]<-NA
    XdifIDCL<-(Xtotal.IDCL-chop.obs.X)
    dif.sq<-(XdifIDCL)^2
    point.pe.IDCL<-sqrt((sum(dif.sq,na.rm=T))/sum(chop.obs.X^2,na.rm=T))
    
    #total.err.IDCL<-sum(Xtotal.IDCL, na.rm=T)-sum(chop.obs.X, na.rm=T)
    total.err.IDCL<-(sum(XdifIDCL, na.rm=T))
    total.err.IDCL<-abs(total.err.IDCL)
    total.pe.IDCL<-total.err.IDCL/sum(chop.obs.X,na.rm=T)
    
    calendar.dif<-get.calendar(XdifIDCL,m,ncut)
    calendar.pe.IDCL<-sqrt((sum(calendar.dif^2,na.rm=T))/sum(chop.obs.X.diag^2,na.rm=T))
    
    if (Plot.box ==TRUE  )
    {  
      ## Boxplots of the cell errors
      
      DCLdif<-get.diagonalcuts(XdifDCL,m,ncut)
      BDCLdif<-get.diagonalcuts(XdifBDCL,m,ncut)
      IDCLdif<-get.diagonalcuts(XdifIDCL,m,ncut)
      Xdif<-cbind(DCLdif,BDCLdif,IDCLdif)
      Xdif[abs(Xdif)==Inf]<-NA
      limitsy<-quantile(Xdif,c(0.05,0.95),na.rm=TRUE)
      boxplot(Xdif,na.rm=T, outline=T,col=c(2,3,4),ylab="prediction minus observed", 
              ylim=limitsy,names=c("DCL","BDCL","IDCL"),
              main=paste(ncut," periods cut",sep=""),cex.main=0.9)
    }
    
    res.point<-data.frame(num.cuts=ncut,point.pe.DCL,point.pe.BDCL,point.pe.IDCL)
    res.calendar<-data.frame(num.cuts=ncut,calendar.pe.DCL,calendar.pe.BDCL,calendar.pe.IDCL)
    res.total<-data.frame(num.cuts=ncut,total.pe.DCL,total.pe.BDCL,total.pe.IDCL)
    if (Tables==TRUE)
    {
      print(round(res.point,num.dec))
      print(round(res.calendar,num.dec))
      print(round(res.total,num.dec))
    }
    pe.vector<-as.numeric(cbind(res.point,res.calendar[,-1],res.total[,-1]))
#     return(list(pe.vector=pe.vector,pe.point=res.point,pe.calendar=res.calendar,
#                 pe.total=res.total,Xdif=Xdif))
    return(list(pe.vector=pe.vector,Xdif=Xdif))
  }
}
