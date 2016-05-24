################################################################################
sim.freq<-function(nbal=4,nbloc=5,nbpop=3,N=1000,mig=0.001,mut=0.0001,f=0){
  #allows for different N and f for each population
  #modified param so that it reflects correctly population effective size
  genofreq<-function(freq,fi){
    if (fi==0) {geno.freq<-outer(freq,freq)}
    else {geno.freq<-outer(freq,freq)*(1-fi)+diag(nbal)*diag(freq)*fi}
    return(geno.freq)
  }
  if (nbal>99) stop ("Too many alleles, must be <100. Exiting")
  cl<-match.call()
  pl<-vector("list",nbloc)
  freq<-rdirichlet(nbloc,rep(1,nbal)) #1=uniform freq has dim nbloc,nbal
  if (length(N)!=nbpop) {
    if (length(N)==1) N<-rep(N,nbpop) 
    else stop("N must be a vector of length nbpop. Exiting.")} 
  if (length(f)!=nbpop){
    if (length(f)==1) f<-rep(f,nbpop) 
    else stop("f must be a vector of length nbpop. Exiting.")}
  param<-outer(4*N/(1+f)*(mig+mut),freq,"*") #verify this [nbpop,nbloc,nbal]
  for (il in 1:nbloc){
    x<-matrix(numeric(nbal*nbpop),nrow=nbpop)
    for (ip in 1:nbpop) x[ip,]<-rdirichlet(1,param[ip,il,])
    pl[[il]]<-x
  }
  gf<-vector("list",nbloc)
  for (il in 1:nbloc){    
    gf[[il]]<-matrix(numeric(nbal^2*nbpop),ncol=nbpop)
    for (ip in 1:nbpop) gf[[il]][,ip]<-genofreq(pl[[il]][ip,],f[ip])}
  if (nbal<10)
    nfun<-function(x,y) x*10+y
  else
    nfun<-function(x,y) x*100+y
  gn<-as.numeric(outer(1:nbal,1:nbal,nfun))
  return(list(call=cl,fpl=pl,gf=gf,gn=gn))
}
#########################################################################
sim.genot<-function(size=50,nbal=4,nbloc=5,nbpop=3,N=1000,mig=0.001,mut=0.0001,f=0){
  a<-sim.freq(nbal,nbloc,nbpop,N,mig,mut,f)
  dat<-data.frame(rep(1:nbpop,each=size),matrix(numeric(nbloc*nbpop*size),ncol=nbloc))
  names(dat)<-c("Pop",paste("loc",1:nbloc,sep="."))
  dumf<-function(x) sample(a$gn,size=size,replace=TRUE,prob=x)
  for (il in 1:nbloc) dat[,il+1]<-as.numeric(apply(a$gf[[il]],2,dumf))
  return(dat)
}
