sim.snp.expsurv.sctest <-
function(n,gtprev,lam,a,b,ztest,diag=FALSE)
  {
    ng=as.integer(rmultinom(1,n,gtprev))
    atime=rexp(n,rep(lam,ng))
    ctime=runif(n,a,b)
    otime=pmin(atime,ctime)
    event=as.integer(atime<ctime)
    SNP=rep(0:2,ng)
    
    SNP[SNP==0]=ztest[1]
    SNP[SNP==1]=ztest[2]
    SNP[SNP==2]=ztest[3]
      
    SNP=factor(SNP)
    
    mod=summary(coxph(Surv(otime,event)~SNP))
    if(diag)
      print(mod)
    c(event=mean(event),pval=mod$sctest[3])
  }
  
  