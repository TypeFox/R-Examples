var.Np<-function(pars,mix.terms,width,z,zdim,pt,n,hessian,N,pa.vec,pa,x){
   # calculate the variance of N and the average p

   # most of this is taken from the summary.ds method in mrds
   # by Jeff Laake, http://github.com/jlaake/mrds

   ret<-list()

   # function to calculate N, to be passed to DeltaMethod()
   Nfct<-function(pars,mix.terms,width,z,zdim,pt,n){
      mu<-mu.calc(pars,mix.terms,width,z,zdim,pt)

      if(length(mu)!=n){
         mu<-rep(mu,n)
      }

      if(!pt){
         return(width*sum(1/mu))
      }else{
         return((pi*width^2)*sum(1/mu))
      }
   }
   # quick function to be passed to DeltaMethod()
   # not really sure how this works, taken from mrds
   pafct<-function(pars,mix.terms,width,z,zdim,pt,n){
      mu<-mu.calc(pars,mix.terms,width,z,zdim,pt)
      if(length(mu)!=n){
         mu<-rep(mu,n)
      }
      return(sum(mu/mu))
   }

   #vcov<-solvecov(hessian)$inv
   vcov<-hessian # this is already inverted in fitmix()

eps<-(.Machine$double.eps)^(1/3)
#######################################################
   fpar<-pars
   fpar1<-pars
   parmat<-NULL
#eps<-.0001

   # Compute first partial (numerically) of log(f(y)) for each observation 
   # for each parameter and store in parmat (n by length(fpar))
   for (i in 1:length(fpar)){
      deltap<-eps*fpar1[i]
      if(deltap==0){
         deltap<-eps
      }
      fpar[i] <- fpar1[i]- deltap
#      x1=-flt.lnl(fpar, ddfobj,TCI,misc.options)
#      x1<- flt(fpar,x,width,mix.terms,0,"hn",z,zdim,pt)
x1<-   log(rowSums(eval.pdf(fpar,x,width,mix.terms,0,"hn",z,zdim,pt)))
      fpar[i] <- fpar1[i]+deltap
#      x2=-flt.lnl(fpar, ddfobj,TCI,misc.options)
#      x2<- flt(fpar,x,width, mix.terms,0,"hn",z,zdim,pt)
x2<-   log(rowSums(eval.pdf(fpar,x,width,mix.terms,0,"hn",z,zdim,pt)))
      parmat=cbind(parmat,(x2-x1)/(2*deltap))
   }

   # Compute varmat using first partial approach (pg 62 of Buckland et al 2002)
   varmat=matrix(0,ncol=length(fpar1),nrow=length(fpar1))

   for(i in 1:length(fpar1)){
      for(j in 1:length(fpar1)){
         varmat[i,j]=sum(parmat[,i]*parmat[,j])
      }
   }

   vcov<-solvecov(varmat)$inv

######################################################

   ### calculate the variance, se and CV of N in the covered region
#   eps<-0.001
#   eps<-sqrt(.Machine$double.eps)

   if(length(pa.vec)==1){
      pa.vec<-rep(pa.vec,n)
   }

   Nhatvar.list<-DeltaMethod(pars,Nfct,vcov,eps,mix.terms=mix.terms,
                             width=width,z=z,zdim=zdim,pt=pt,n=n)

   Nhatvar<-Nhatvar.list$variance + sum((1-pa.vec)/pa.vec^2)
   # put them in the results list
   ret$N.var<-Nhatvar
   ret$N.se<-sqrt(Nhatvar)
   ret$N.CV<-sqrt(Nhatvar)/N
   cvN<-sqrt(Nhatvar)/N

   ### calculate variance, se and CV for average detectability

   vc1.list<-DeltaMethod(pars,pafct,vcov,eps,mix.terms=mix.terms,
                      width=width,z=z,zdim=zdim,pt=pt,n=n)
   vc1<-vc1.list$variance
   vc2<-sum(1 - pa.vec)
   covar<-sum((1 - pa.vec)/pa.vec)
   var.pbar.list<-list(var=vc1+vc2,partial=vc1.list$partial,covar=covar)

   covar<-t(Nhatvar.list$partial)%*%vcov%*%var.pbar.list$partial+
                        var.pbar.list$covar
   var.pbar<-pa^2*(cvN^2 + var.pbar.list$var/n^2-2*covar/(n*N))
   # put the results in the return list...
   ret$pa.var<-var.pbar
   ret$pa.se<-sqrt(var.pbar)
   ret$pa.CV<-sqrt(var.pbar)/pa
   ret$vcov<-vcov

   return(ret)
}
