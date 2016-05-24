hap.control <- function(mb=0,pr=0,po=0.001,to=0.001,th=1,maxit=100,n=0,
      ss=0,rs=0,rp=0,ro=0,rv=0,sd=0,mm=0,mi=0,mc=50,ds=0.1,de=0,q=0,
      hapfile="hap.out",assignfile="assign.out")
{
   list(mb=mb,pr=pr,po=po,to=to,th=th,maxit=maxit,n=n,ss=ss,rs=rs,
       rp=rp,ro=ro,rv=rv,sd=sd,mm=mm,mi=mi,mc=mc,ds=ds,de=de,q=q,
       hapfile=hapfile,assignfile=assignfile)
}

hap<-function(id,data,nloci,loci=rep(2,nloci),names=paste("loci",1:nloci,sep=""),
              control=hap.control())
{
  if (control$rv & control$ro) stop("rv and ro flags cannot both be set in hap.control\n");
# if (mi==0 & (mc | ds | de)) stop("mc, ds, de parameters are only legal if mi is set\n");
  if (control$rp & control$mm==0) stop("rp option only relevant with mm # option\n");
  nobs<-dim(data)[1]
  data<-as.matrix(data)
  if(length(id)!=dim(data)[1]) stop("id and data should have the same length")
  l1<-niter<-converge<-0

  z<-.C("hap",nobs=as.integer(nobs),idstr=as.character(id),data=as.character(t(data)),
        nloci=as.integer(nloci),loci=as.integer(loci),names=as.character(names),
        mb=as.double(control$mb),
        pr=as.double(control$pr),
        po=as.double(control$po),
        to=as.double(control$to),
        th=as.double(control$th),
        maxitt=as.double(control$maxit),
        n=as.integer(control$n),
        sst=as.integer(control$ss),
        rst=as.integer(control$rs),
        rp=as.integer(control$rp),
        ro=as.integer(control$ro),
        rv=as.integer(control$rv),
        sd=as.double(control$sd),
        mm=as.integer(control$mm),
        mi=as.integer(control$mi),
        mc=as.integer(control$mc),
        ds=as.double(control$ds),
        de=as.double(control$de),
        q=as.integer(control$q),
        l1=as.double(l1),
        niter=as.integer(niter),
        converged=as.integer(converge),
        hapfile=as.character(control$hapfile),
        assignfile=as.character(control$assignfile),
        PACKAGE="gap")

  list(l1=z$l1,converge=z$converged,niter=z$niter)
}
