ce.dat.prep <-
function(xt.dat, failure.dat, ref_time=NULL){ 
  # modify the case when the failure time and the maximum covariate time cannot match
  # very slow 
  names(xt.dat)[1:2]=c("id", "time")
  
  max.t.obj = aggregate(xt.dat$time, by=list(id=xt.dat$id), max)[failure.dat[,"status"]==1,]
  nomatch.id = max.t.obj$id[max.t.obj$x<failure.dat[failure.dat[,"status"]==1, "time"]]
  nomatch.dat.id = match(interaction(max.t.obj[nomatch.id,]),interaction(xt.dat[,1:2]))
  
  dat=rbind(xt.dat, cbind(nomatch.id, failure.dat[nomatch.id, "time"], 
                          xt.dat[nomatch.dat.id,-c(1,2)]))
  
  
  npts.vec=as.numeric(table(as.factor(xt.dat$id)))
  npts=max(npts.vec)
  n=length(unique(xt.dat$id))
  
  p=ncol(xt.dat)-2 # number of dynamic covariates
  
  x.val=matrix(0, n, npts*p)
  wts.mat=matrix(0, n, npts)
  
  if(is.null(ref_time)){
    ref_time=rep(0,n)
  }
  
  for(i in 1:n){
    for(j in 1:p){
      a=(j-1)*npts+1
      b=(j-1)*npts+npts.vec[i]
      x.val[i, a:b]=xt.dat[xt.dat$id==i,j+2]    
    }
    
    wts.mat[i, 1:npts.vec[i]]=diff(c(ref_time[i], xt.dat$time[xt.dat$id==i]))
  }
  
  failure.dat=as.data.frame(cbind(failure.dat[,"time"], failure.dat[,"status"]))
  names(failure.dat)=c("failure.time","delta")
  
  aux.inf=list(wts.mat=wts.mat, npts.vec=npts.vec)
  xt.obj=list(x.val=x.val, npts=npts, n=n)
  ce.obj=list(failure.dat=failure.dat, xt.obj=xt.obj, aux.inf=aux.inf)
  
  class(ce.obj)="ce.obj"
  return(ce.obj)
}
