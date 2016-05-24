meta.RiskR <-
function(data.mi, BB.grdnum=1000, B.sim=20000, cov.prob=0.95, print=T, studyCI=T, midp=T, ratio.upper=1000)
  {
    data.original=data.mi
  
    if(sum(data.mi[,1])+sum(data.mi[,2])>0){    
    data.mi=data.mi[data.mi[,1]+data.mi[,2]>0, ,drop=F]
    n=length(data.mi[,1])
    d1=data.mi[,1];             d2=data.mi[,2]
    n1=data.mi[,3];             n2=data.mi[,4]
    p1=d1/n1;                   p2=d2/n2
    id=(1:n)[d1*d2==0]
    p1[id]=(d1[id]+0.5)/(n1[id]+1);  p2[id]=(d2[id]+0.5)/(n2[id]+1)

    weight=(n1*n2/(n1+n2))/sum(n1*n2/(n1+n2))
    
    varp1=p1*(1-p1)/n1
    varp2=p2*(1-p2)/n2
    p1.pool=sum(p1*weight)
    p2.pool=sum(p2*weight)

    mu.MH=log(p2.pool)-log(p1.pool)
    sd.MH=sqrt(sum(varp2*weight^2)/p2.pool^2+sum(varp1*weight^2)/p1.pool^2)    
    ci.MH=c(mu.MH-qnorm((1+cov.prob)/2)*sd.MH, mu.MH+qnorm((1+cov.prob)/2)*sd.MH)
    p.MH=1-pchisq(mu.MH^2/sd.MH^2,1)
 
    d0=max(abs(ci.MH))
    BB.grdnum=2*round(BB.grdnum/2)+1
    delta.grd=exp(seq(-min(log(ratio.upper), d0*4), min(log(ratio.upper), d0*4),length=BB.grdnum))
    

    pv1.pool=pv2.pool=numeric(0)
    for(kk in 1:n)
      { xx1=data.mi[kk,1]
        xx2=data.mi[kk,2]
        nn1=data.mi[kk,3] 
        nn2=data.mi[kk,4]
        fit=priskR.exact(xx1,xx2,nn1,nn2, delta.grd)
        pv1.pool=rbind(pv1.pool, fit$pv1); pv2.pool=rbind(pv2.pool, fit$pv2)
        if(print==T)  cat("study=", kk, "\n")
      }


    for(i in 1:n)
    for(j in 1:BB.grdnum)
        pv1.pool[i,(BB.grdnum-j+1)]=max(pv1.pool[i,1:(BB.grdnum-j+1)]);pv2.pool[i,j]=max(pv2.pool[i,j:BB.grdnum])
          
    sigma0=1/n1+1/n2
   
    set.seed(100)
    tnull1=tnull2=matrix(0,B.sim,3)
    y=matrix(runif(B.sim*n), n, B.sim)
    y=y/(1+1e-2)
    tnull1[,1]=apply(-log(1-y)/sigma0*(d2>0), 2, sum)
    tnull1[,2]=apply(y/sigma0*(d2>0), 2, sum)
    tnull1[,3]=apply(asin(y)/sigma0*(d2>0), 2, sum)
    tnull2[,1]=apply(-log(1-y)/sigma0*(d1>0), 2, sum)
    tnull2[,2]=apply(y/sigma0*(d1>0), 2, sum)
    tnull2[,3]=apply(asin(y)/sigma0*(d1>0), 2, sum)
       
    if(sum(d1)*sum(d2)!=0)  alpha0=(1+cov.prob)/2;
    if(sum(d1)*sum(d2)==0)  alpha0=cov.prob
 
    cut1=cut2=rep(0,3)
    for(b in 1:3)
       {cut1[b]=quantile(tnull1[,b], 1-alpha0)
        cut2[b]=quantile(tnull2[,b], 1-alpha0)
        }

    t1=t2=matrix(0,BB.grdnum,3)
    pv1.pool=pv1.pool/(1+1e-2)
    pv2.pool=pv2.pool/(1+1e-2)
    t1[,1]=apply(-log(1-pv1.pool)/sigma0*(d2>0), 2, sum);  t2[,1]=apply(-log(1-pv2.pool)/sigma0*(d1>0), 2, sum)
    t1[,2]=apply(pv1.pool/sigma0*(d2>0), 2, sum);          t2[,2]=apply(pv2.pool/sigma0*(d1>0), 2, sum)
    t1[,3]=apply(asin(pv1.pool)/sigma0*(d2>0), 2, sum);    t2[,3]=apply(asin(pv2.pool)/sigma0*(d1>0), 2, sum)
    

    ci.fisher=  c(min(delta.grd[t1[,1]>=cut1[1]]),max(delta.grd[t2[,1]>=cut2[1]]))
    ci.cons=    c(min(delta.grd[t1[,2]>=cut1[2]]),max(delta.grd[t2[,2]>=cut2[2]]))
    ci.iv=c(min(delta.grd[t1[,3]>=cut1[3]]),max(delta.grd[t2[,3]>=cut2[3]]))    
    ci.MH=exp(ci.MH)
    ci.range=c(min(delta.grd), max(delta.grd))


    est.fisher=delta.grd[abs(t2[,1]-t1[,1])==min(abs(t2[,1]-t1[,1]))][1]    
    est.cons=delta.grd[abs(t2[,2]-t1[,2])==min(abs(t2[,2]-t1[,2]))][1]    
    est.iv=delta.grd[abs(t2[,3]-t1[,3])==min(abs(t2[,3]-t1[,3]))][1]    
    est.MH=exp(mu.MH)
    est.range=NA

    if(sum(d1)==0)
      {ci.fisher[2]=ci.cons[2]=ci.iv[2]=Inf;
       est.fisher=est.cons=est.iv=Inf}

    if(sum(d2)==0)
      {ci.fisher[1]=ci.cons[1]=ci.iv[1]=0
       est.fisher=est.cons=est.iv=0}

    n0=(BB.grdnum+1)/2
    c1=t1[n0,]; c2=t2[n0,]

    p.fisher=  min(1, 2*min(c(1-mean(tnull1[,1]>=c1[1]), 1-mean(tnull2[,1]>=c2[1]))))
    p.cons=    min(1, 2*min(c(1-mean(tnull1[,2]>=c1[2]), 1-mean(tnull2[,2]>=c2[2]))))
    p.iv=min(1, 2*min(c(1-mean(tnull1[,3]>=c1[3]), 1-mean(tnull2[,3]>=c2[3]))))


    pvalue=c(p.cons, p.iv, p.fisher, p.MH, NA)
    ci=cbind(ci.cons, ci.iv, ci.fisher,ci.MH, ci.range)
    ci=rbind(c(est.cons, est.iv, est.fisher, est.MH, est.range), ci, pvalue)
    rownames(ci)=c("est", "lower CI", "upper CI", "p")
    colnames(ci)=c("constant", "inverse-variance", "fisher", "asymptotical-MH", "range")
    }else
    {ci=rbind(rep(NA, 4), rep(0, 4), rep(Inf, 4), rep(1, 4))
     rownames(ci)=c("est", "lower CI", "upper CI", "p")
     colnames(ci)=c("constant", "inverse-variance", "fisher", "asymptotical-MH")
    }

  
###################################################################################################  

    study.ci=NULL

    if(studyCI==T)
      {data.mi=data.original

       n=length(data.mi[,1])
       study.ci=matrix(0, n, 5)
       colnames(study.ci)=c("est", "lower CI", "upper CI", "p",  "limit")
       rownames(study.ci)=1:n

       for(kk in 1:n) 
          {xx1=data.mi[kk,1]
           xx2=data.mi[kk,2] 
           nn1=data.mi[kk,3] 
           nn2=data.mi[kk,4]

           fit=ci.RiskR(xx1, xx2, nn1, nn2, cov.prob=cov.prob, BB.grdnum=BB.grdnum, midp=midp)
           study.ci[kk,2]=fit$lower
           study.ci[kk,3]=fit$upper
           study.ci[kk,1]=fit$est
           study.ci[kk,4]=fit$p
           study.ci[kk,5]=fit$status
           
           rownames(study.ci)[kk]=paste("study ", kk)
           }
       }       
 

    return(list(ci.fixed=ci, study.ci=study.ci, precision=paste("+/-", (max(log(delta.grd))-min(log(delta.grd)))/BB.grdnum)))
 }
