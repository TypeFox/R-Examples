meta.RateD <-
function(data.mi, BB.grdnum=1000, B.sim=10000, cov.prob=0.95, midp=T, MH.imputation=F, print=T, studyCI=T)
  {
    n=length(data.mi[,1])
    e1=data.mi[,3];             e2=data.mi[,4]
    lambda1=data.mi[,1]/data.mi[,3]; lambda2=data.mi[,2]/data.mi[,4]
 
    if(MH.imputation==T)
      {id=(1:n)[lambda1*lambda2==0]
       lambda1[id]=(data.mi[id,1]+0.5)/e1[id];  lambda2[id]=(data.mi[id,2]+0.5)/e2[id]
       }

    deltalambda=lambda2-lambda1
    varlambda=lambda2/e2+lambda1/e1
    weight=(e1*e2/(e1+e2))/sum(e1*e2/(e1+e2))
    mu.MH=sum(deltalambda*weight);  sd.MH=sqrt(sum(weight^2*varlambda))
    ci.MH=c(mu.MH-qnorm((1+cov.prob)/2)*sd.MH, mu.MH+qnorm((1+cov.prob)/2)*sd.MH)
    p.MH=1-pchisq(mu.MH^2/sd.MH^2,1)

 
    d0=max(abs(ci.MH))

    BB.grdnum=2*(round(BB.grdnum/2))+1
    delta.grd=seq(-d0*5, d0*5,length=BB.grdnum); 
    
    pv1.pool=pv2.pool=numeric(0)
    for(kk in 1:n)
      { xx1=data.mi[kk,1]
        xx2=data.mi[kk,2] 
        ee1=data.mi[kk,3] 
        ee2=data.mi[kk,4]
        fit=prateD.exact(xx1, xx2, ee1, ee2, delta.grd, midp=midp)
        pv1.pool=rbind(pv1.pool, fit$pv1); pv2.pool=rbind(pv2.pool, fit$pv2)
        if(print==T)
          cat("study=", kk, "\n")
      }


    for(i in 1:n)
      { for(j in 1:BB.grdnum)
          { pv1.pool[i,(BB.grdnum-j+1)]=max(pv1.pool[i,1:(BB.grdnum-j+1)]);pv2.pool[i,j]=max(pv2.pool[i,j:BB.grdnum])
          }
      }
  

    sigma0=1/data.mi[,3]+1/data.mi[,4]
   
    set.seed(100)
    tnull=matrix(0,B.sim,3)
    y=matrix(runif(B.sim*n), n, B.sim)
    y=y/(1+1e-2)
    tnull[,1]=apply(-log(1-y)/sigma0, 2, sum)
    tnull[,2]=apply(y/sigma0, 2, sum)
    tnull[,3]=apply(asin(y)/sigma0, 2, sum)


    alpha0=(1+cov.prob)/2; 
    cut=rep(0,3)
    for(b in 1:3)
       cut[b]=quantile(tnull[,b], 1-alpha0)
        

    t1=t2=matrix(0,BB.grdnum,3)
    pv1.pool=pv1.pool/(1+1e-2)
    pv2.pool=pv2.pool/(1+1e-2)
    t1[,1]=apply(-log(1-pv1.pool)/sigma0, 2, sum);  t2[,1]=apply(-log(1-pv2.pool)/sigma0, 2, sum)
    t1[,2]=apply(pv1.pool/sigma0, 2, sum);          t2[,2]=apply(pv2.pool/sigma0, 2, sum)
    t1[,3]=apply(asin(pv1.pool)/sigma0, 2, sum);    t2[,3]=apply(asin(pv2.pool)/sigma0, 2, sum)
    

    ci.fisher=  c(min(delta.grd[t1[,1]>=cut[1]]),max(delta.grd[t2[,1]>=cut[1]]))
    ci.cons=    c(min(delta.grd[t1[,2]>=cut[2]]),max(delta.grd[t2[,2]>=cut[2]]))
    ci.iv=c(min(delta.grd[t1[,3]>=cut[3]]),max(delta.grd[t2[,3]>=cut[3]]))    
    ci.MH=ci.MH
    ci.range=c(min(delta.grd), max(delta.grd))

    est.fisher=delta.grd[abs(t2[,1]-t1[,1])==min(abs(t2[,1]-t1[,1]))][1]    
    est.cons=delta.grd[abs(t2[,2]-t1[,2])==min(abs(t2[,2]-t1[,2]))][1]    
    est.iv=delta.grd[abs(t2[,3]-t1[,3])==min(abs(t2[,3]-t1[,3]))][1]    
    est.MH=exp(mu.MH)
    est.range=NA


    n0=(BB.grdnum+1)/2
    c1=t1[n0,]; c2=t2[n0,]

    p.fisher=  min(1, 2*min(c(1-mean(tnull[,1]>=c1[1]), 1-mean(tnull[,1]>=c2[1]))))
    p.cons=    min(1, 2*min(c(1-mean(tnull[,2]>=c1[2]), 1-mean(tnull[,2]>=c2[2]))))
    p.iv=min(1, 2*min(c(1-mean(tnull[,3]>=c1[3]), 1-mean(tnull[,3]>=c2[3]))))


    pvalue=c(p.cons, p.iv, p.fisher, p.MH, NA)
    ci=cbind(ci.cons, ci.iv, ci.fisher,ci.MH, ci.range)
    ci=rbind(c(est.cons, est.iv, est.fisher, est.MH,  est.range), ci, pvalue)
    rownames(ci)=c("est", "lower CI", "upper CI", "p")
    colnames(ci)=c("constant", "inverse-variance", "fisher", "asymptotical-MH", "range")
  
###########################################################################################
  
   study.ci=NULL
      
   if(studyCI==T)
     {
      n=length(data.mi[,1])
      study.ci=matrix(0, n, 5)
      colnames(study.ci)=c("est", "lower CI", "upper CI", "p", "limit")
      rownames(study.ci)=1:n
      for(kk in 1:n) 
         {xx1=data.mi[kk,1]
          xx2=data.mi[kk,2] 
          ee1=data.mi[kk,3] 
          ee2=data.mi[kk,4]

          fit=ci.RateD(xx1, xx2, ee1, ee2, cov.prob=cov.prob, BB.grdnum=BB.grdnum, midp=midp)
          study.ci[kk,2]=fit$lower
          study.ci[kk,3]=fit$upper
          study.ci[kk,1]=fit$est
          study.ci[kk,4]=fit$p
          study.ci[kk,5]=fit$status

          rownames(study.ci)[kk]=paste("study ", kk)
          }
       }


    return(list(ci.fixed=ci, study.ci=study.ci, precision=paste("+/-", (max(delta.grd)-min(delta.grd))/BB.grdnum)))

 }
