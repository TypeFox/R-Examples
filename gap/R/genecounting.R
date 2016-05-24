# 10/4/2005, 13/4/2005
gc.control <- function(xdata=FALSE, convll=1,handle.miss=0,eps=0.000001,
                       tol=0.00000001, maxit=50,pl=0.001,assignment="assign.dat",verbose=T)
{
   list(xdata=xdata,convll=convll,handle.miss=handle.miss,eps=eps,tol=tol,
        maxit=maxit,pl=pl,assignment=assignment,verbose=verbose)
}

genecounting <- function(data,weight=NULL,loci=NULL,control=gc.control())
{
  if (control$xdata)
  {
     sex <- data[,1]
     data <- data[,-1]
  }
  else
  sex <- rep(2,dim(data)[1])
  if(is.null(weight)) weight<-rep(1,dim(data)[1])
# precis<-1
# to call dpmach
# tol<-1.2
#  while(tol>1.0)
#  {
#     precis<-precis/2.0;
#     tol<-1.0+precis;
#  }
  precis<-.Machine$double.eps
  gid<-1:(dim(data)[1])
  nloci=dim(data)[2]/2
  if(is.null(loci))
  {
    loci<-rep(0,nloci)
    for (i in 1:nloci)
    {
        loci[i]=max(data[,c(2*i-1,2*i)],na.rm=TRUE)
    }
  }
  data<-as.matrix(data)
  data<-t(data)
  hapall<-1
  for(i in 1:nloci)
  {
    hapall<-hapall*loci[i];
  }
  h0<-h1<-hapid<-rep(0,hapall)
  obscom<-length(weight)
  prob<-rep(0,obscom)
  lnl0<-lnl1<-0
  npusr<-npdat<-rep(0,2)
# 13/11/2003
# change to reduce memory request
# htrtable<-matrix(rep(0,obscom*hapall),nrow=obscom)
  verbose<-0
  iter<-0
  converge<-0
  z <- .C("gcx",verbose=as.integer(control$verbose),
           Rhandlemissing=as.integer(control$handle.miss),
           Rconvll=as.integer(control$convll),
           Reps=as.double(control$eps),
           Rtol=as.double(control$tol),
           Rmaxit=as.integer(control$maxit),
           Rpl=as.double(control$pl),
           precis=as.double(precis),
           gid=as.integer(gid),
           Rnloci=as.integer(nloci),
           Rloci=as.integer(loci),
           Robscom=as.integer(obscom),
           Rhapall=as.integer(hapall),
           genotype=as.integer(data),
           count=as.integer(weight),
           Rxdata=as.integer(control$xdata),
           sex=as.integer(sex),
           hapid=as.integer(hapid),
           prob=as.double(prob),
           Rh0=as.double(h0),
           Rh1=as.double(h1),
           lnl0=as.double(lnl0),
           lnl1=as.double(lnl1),
           npusr=as.integer(npusr),
           npdat=as.integer(npdat),
           iter=as.integer(iter),
           converge=as.integer(converge),assignment=as.character(control$assignment),PACKAGE="gap"
           )
  x<-0
  hapid<-0
# Dprime<-sum(z$Rh0*abs(z$Rh1-z$Rh0))
# x<-t(matrix(z$htrtable/2,nrow=hapall))
# hapid<-apply(x,2,sum)>0
# x<-x[,hapid]
# hapid<-(1:hapall)[hapid]
  di0<-1-sum((z$Rh0)^2)
  di1<-1-sum((z$Rh1)^2)
  prob<-z$prob/sum(z$prob)
  resid<-weight*(1-prob)
  list(h=z$Rh1, h0=z$Rh0, prob=prob, l0=z$lnl0, l1=z$lnl1,
       hapid=hapid, npusr=z$npusr, npdat=z$npdat, htrtable=x,
       iter=z$iter,converge=z$converge,di0=di0,di1=di1,resid=resid)
}
