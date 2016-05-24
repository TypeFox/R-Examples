backfit.lm.main=function(x, y, cutp, mincut=0)
{n=length(y)
 p=length(x[1,])
 

 if(mincut>0)
   {effect.range=ceiling(n*mincut):floor(n*(1-mincut))}
 if(mincut==0)
   {effect.range=1:n}


 order.index=matrix(0, n, p)
 for(i in 1:p)
    {order.index[,i]=order(x[,i])}


 muhat=mean(y)
 sigma=sum((y-muhat)^2)/(n-1)

 y.pool=x.replicate=matrix(0, n, p)
 for(i in 1:p)
    {y.pool[,i]=y[order.index[,i]]
     x.replicate[,i]=rep((1:n)[cumsum(table(x[,i]))], table(x[,i]))}

 
 score.mat=t((t(x)<cutp)*1)
 score.tot=apply(score.mat,1,sum)       

 zscore=sum((y-muhat)*score.tot)/sqrt(sigma*(sum(score.tot^2)-(sum(score.tot))^2/n))
     
 count=1
 loop=1
 while(loop==1)
      {
       score.mat=t((t(x)<cutp)*1)
       score.tot=apply(score.mat,1,sum)       
   
       score.pool=muhat.pool=matrix(0, n, p)
       for(i in 1:p)
          {idx=order.index[,i] 
           score.pool[,i]=(score.tot-score.mat[,i])[idx]
           muhat.pool[,i]=muhat[idx]}
       
       score.stat=t(t(apply(y.pool-muhat,2,cumsum))+apply((y.pool-muhat)*score.pool,2,sum))
       v.stat=sigma*(t(t(2*apply(score.pool,2,cumsum)+1:n)+apply(score.pool^2,2,sum))-(t(matrix(rep(apply(score.pool,2,sum),n),p,n))+1:n)^2/n)

       for(i in 1:p)
          {idx=x.replicate[,i]
           v.stat[,i]=v.stat[idx ,i]   
           score.stat[,i]=score.stat[idx, i]
          } 

       test=score.stat/sqrt(v.stat+1e-8)

       if(mincut==0)
         {mtest=apply(test,2,max)}

       if(mincut>0)
         {mtest=rep(0, p)
          for(i in 1:p)
             {
              test.post=test[max(effect.range)+1,i]
              test0=test[effect.range,i]
              test0=test0[test0!=test.post]
              mtest[i]=max(test0)}
          }
       loop=0
       if(max(mtest)>zscore[count]+1e-8)
         {loop=1
          i=(1:p)[mtest==max(mtest)][1]
          i0=max(effect.range[test[effect.range,i]==mtest[i]])+1
          if(i0>n)
            mcut=Inf
          if(i0<=n)
            mcut=x[order.index[,i],i][i0]
          cutp[i]=mcut
          zscore=c(zscore, max(mtest)) 
          count=count+1           
         }
        
       }
       
       
  return(list(cutp=cutp, zscore=zscore))
 }

############################################################################################################################################
############################################################################################################################################
############################################################################################################################################

backfit.lm.interaction=function(x, trt, y, cutp, mincut=0)
{n=length(y)
 p=length(x[1,])
 

 if(mincut>0)
   {effect.range=ceiling(n*mincut):floor(n*(1-mincut))}
 if(mincut==0)
   {effect.range=1:n}


 order.index=matrix(0, n, p)
 for(i in 1:p)
    {order.index[,i]=order(x[,i])}

 fit=lm(y~trt)
 muhat=fit$fitt
 sigma=sum((y-muhat)^2)/(n-2)
 sigma.sum=sum(solve(t(cbind(1,trt))%*%cbind(1,trt))*sigma)

 y.pool=trt.pool=muhat.pool=x.replicate=matrix(0, n, p)
 for(i in 1:p)
    {y.pool[,i]=y[order.index[,i]]
     trt.pool[,i]=trt[order.index[,i]]
     muhat.pool[,i]=muhat[order.index[,i]]
     x.replicate[,i]=rep((1:n)[cumsum(table(x[,i]))], table(x[,i]))}

 
 score.mat=t((t(x)<cutp)*1)
 score.tot=apply(score.mat,1,sum)       
 zscore=sum((y-muhat)*score.tot*trt)/sqrt(sigma*sum(score.tot^2*trt)-(sum(score.tot*trt))^2*sigma.sum)
     
 count=1
 loop=1
 while(loop==1)
      {
       score.mat=t((t(x)<cutp)*1)
       score.tot=apply(score.mat,1,sum)       
   
       score.pool=matrix(0, n, p)
       for(i in 1:p)
          {idx=order.index[,i] 
           score.pool[,i]=(score.tot-score.mat[,i])[idx]}

       score.stat=t(t(apply((y.pool-muhat.pool)*trt.pool,2,cumsum))+apply((y.pool-muhat.pool)*score.pool*trt.pool,2,sum))
       v11.stat=sigma*t(t(apply((2*score.pool+1)*trt.pool,2,cumsum))+apply(score.pool^2*trt.pool,2,sum))
       v22.stat=t((apply(score.pool*trt.pool,2,sum)+t(apply(trt.pool,2,cumsum)))^2)*sigma.sum
       v.stat=(v11.stat-v22.stat)


       for(i in 1:p)
          {idx=x.replicate[,i]
           v.stat[,i]=v.stat[idx ,i]   
           score.stat[,i]=score.stat[idx, i]
          } 

       test=score.stat/sqrt(v.stat+1e-8)

       if(mincut==0)
         {mtest=apply(test,2,max)}

       if(mincut>0)
         {mtest=rep(0, p)
          for(i in 1:p)
             {
              test.post=test[max(effect.range)+1,i]
              test0=test[effect.range,i]
              test0=test0[test0!=test.post]
              mtest[i]=max(test0)}
          }
       loop=0
       if(max(mtest)>zscore[count]+1e-8)
         {loop=1
          i=(1:p)[mtest==max(mtest)][1]
          i0=max(effect.range[test[effect.range,i]==mtest[i]])+1
          if(i0>n)
            mcut=Inf
          if(i0<=n)
            mcut=x[order.index[,i],i][i0]
          cutp[i]=mcut
          zscore=c(zscore, max(mtest)) 
          count=count+1           
         }
        
       }
       
       
  return(list(cutp=cutp, zscore=zscore))
 }

############################################################################################################################################
############################################################################################################################################
############################################################################################################################################

backfit.logistic.main=function(x, y, cutp, mincut=0, wt=wt)
{n=length(y)
 p=length(x[1,])
 

 if(mincut>0)
   {effect.range=ceiling(n*mincut):floor(n*(1-mincut))}
 if(mincut==0)
   {effect.range=1:n}


 order.index=matrix(0, n, p)
 for(i in 1:p)
    {order.index[,i]=order(x[,i])}

 wt.pool=y.pool=x.replicate=matrix(0, n, p)
 for(i in 1:p)
    {y.pool[,i]=y[order.index[,i]]
     wt.pool[,i]=wt[order.index[,i]]
     x.replicate[,i]=rep((1:n)[cumsum(table(x[,i]))], table(x[,i]))}

 phat=sum(y*wt)/sum(wt)
 sigma=as.numeric(summary(glm(y~1, family="binomial", weight=wt))$cov.un*(phat*(1-phat))^2*n)

 score.mat=t((t(x)<cutp)*1)
 score.tot=apply(score.mat,1,sum)       

 zscore=sum((y-phat)*score.tot)/sqrt(phat*(1-phat)*(sum(score.tot^2)-(sum(score.tot))^2/n))
     
 count=1
 loop=1
 while(loop==1)
      {
       score.mat=t((t(x)<cutp)*1)
       score.tot=apply(score.mat,1,sum)       
   
       score.pool=matrix(0, n, p)
       for(i in 1:p)
          {idx=order.index[,i] 
           score.pool[,i]=(score.tot-score.mat[,i])[idx]}

       score.stat=t(t(apply((y.pool-phat)*wt.pool,2,cumsum))+apply((y.pool-phat)*(score.pool*wt.pool),2,sum))
       v.stat=sigma*(t(t(2*apply(score.pool*wt.pool^2,2,cumsum)+apply(wt.pool^2,2,cumsum))+apply(score.pool^2*wt.pool^2,2,sum))-(t(matrix(rep(apply(score.pool*wt.pool,2,sum),n),p,n))+apply(wt.pool,2,cumsum))^2/n)

       for(i in 1:p)
          {idx=x.replicate[,i]
           v.stat[,i]=v.stat[idx ,i]   
           score.stat[,i]=score.stat[idx, i]
          } 

       test=score.stat/sqrt(v.stat+1e-8)

       if(mincut==0)
         {mtest=apply(test,2,max)}

       if(mincut>0)
         {mtest=rep(0, p)
          for(i in 1:p)
             {
              test.post=test[max(effect.range)+1,i]
              test0=test[effect.range,i]
              test0=test0[test0!=test.post]
              mtest[i]=max(test0)}
          }
       loop=0
       if(max(mtest)>zscore[count]+1e-8)
         {loop=1
          i=(1:p)[mtest==max(mtest)][1]
          i0=max(effect.range[test[effect.range,i]==mtest[i]])+1
          if(i0>n)
            mcut=Inf
          if(i0<=n)
            mcut=x[order.index[,i],i][i0]
          cutp[i]=mcut
          zscore=c(zscore, max(mtest)) 
          count=count+1           
         }
        
       }
       
       
  return(list(cutp=cutp, zscore=zscore))
 }

############################################################################################################################################
############################################################################################################################################
############################################################################################################################################

backfit.logistic.interaction=function(x, trt, y, cutp, mincut=0, wt=wt)
{n=length(y)
 p=length(x[1,])
 

 if(mincut>0)
   {effect.range=ceiling(n*mincut):floor(n*(1-mincut))}
 if(mincut==0)
   {effect.range=1:n}


 order.index=matrix(0, n, p)
 for(i in 1:p)
    {order.index[,i]=order(x[,i])}

 wt.pool=trt.pool=y.pool=x.replicate=matrix(0, n, p)
 for(i in 1:p)
    {y.pool[,i]=y[order.index[,i]]
     wt.pool[,i]=wt[order.index[,i]]
     trt.pool[,i]=trt[order.index[,i]]
     x.replicate[,i]=rep((1:n)[cumsum(table(x[,i]))], table(x[,i]))}

 fit=glm(y~trt, family="binomial", weight=wt)
 phat=fit$fitt
 sigma=sum(summary(fit)$cov.scale)

 score.mat=t((t(x)<cutp)*1)
 score.tot=apply(score.mat,1,sum)       

 zscore=sum((y-phat)*trt*score.tot)/sqrt(sum(score.tot^2*trt*phat*(1-phat))-(sum(score.tot*trt*phat*(1-phat)))^2*sigma)
     
 count=1
 loop=1
 while(loop==1)
      {score.mat=t((t(x)<cutp)*1)
       score.tot=apply(score.mat,1,sum)       
   
       score.pool=phat.pool=matrix(0, n, p)
       for(i in 1:p)
          {idx=order.index[,i] 
           score.pool[,i]=(score.tot-score.mat[,i])[idx]
           phat.pool[,i]=phat[idx]}

       score.stat=t(t(apply(trt.pool*wt.pool*(y.pool-phat.pool),2,cumsum))+apply((y.pool-phat.pool)*trt.pool*score.pool*wt.pool,2,sum))
       v11.stat=t(t(2*apply(score.pool*trt.pool*wt.pool^2*phat.pool*(1-phat.pool),2,cumsum)+apply(trt.pool*wt.pool^2*phat.pool*(1-phat.pool),2,cumsum))+apply(score.pool^2*wt.pool^2*trt.pool*phat.pool*(1-phat.pool),2,sum))
       v22.stat=(t(t(apply(trt.pool*wt.pool*phat.pool*(1-phat.pool),2,cumsum))+apply(score.pool*trt.pool*wt.pool*phat.pool*(1-phat.pool),2,sum)))^2*sigma
       v.stat=(v11.stat-v22.stat)

       for(i in 1:p)
          {idx=x.replicate[,i]
           v.stat[,i]=v.stat[idx ,i]   
           score.stat[,i]=score.stat[idx, i]
          } 

       test=score.stat/sqrt(v.stat+1e-8)

       if(mincut==0)
         {mtest=apply(test,2,max)}

       if(mincut>0)
         {mtest=rep(0, p)
          for(i in 1:p)
             {
              test.post=test[max(effect.range)+1,i]
              test0=test[effect.range,i]
              test0=test0[test0!=test.post]
              mtest[i]=max(test0)}
          }
       loop=0
       if(max(mtest)>zscore[count]+1e-8)
         {loop=1
          i=(1:p)[mtest==max(mtest)][1]
          i0=max(effect.range[test[effect.range,i]==mtest[i]])+1
          if(i0>n)
            mcut=Inf
          if(i0<=n)
            mcut=x[order.index[,i],i][i0]
          cutp[i]=mcut
          zscore=c(zscore, max(mtest)) 
          count=count+1           
         }
        
       }
       
       
  return(list(cutp=cutp, zscore=zscore))
 }



############################################################################################################################################
############################################################################################################################################
############################################################################################################################################

backfit.cox.main=function(x, y, delta, cutp, mincut=0)
{n=length(y)
 p=length(x[1,])
 id=order(y)
 y=y[id]
 delta=delta[id]
 x=x[id,]


 if(mincut>0)
   {effect.range=ceiling(n*mincut):floor(n*(1-mincut))}
 if(mincut==0)
   {effect.range=1:n}


 order.index=matrix(0, n, p)
 for(i in 1:p)
    {order.index[,i]=order(x[,i])}

 delta.pool=x.replicate=matrix(0, n, p)
 for(i in 1:p)
    {delta.pool[,i]=delta[order.index[,i]]
     x.replicate[,i]=rep((1:n)[cumsum(table(x[,i]))], table(x[,i]))}

 
 ntime=length(unique(y))
 time.index=rep(c(1, cumsum(table(y))+1)[-(ntime+1)], table(y))
 num.risk=(n:1)[time.index]
 time.index.plus=rep(cumsum(table(y)), table(y))

 risk.set=matrix(0, n, n)
 for(i in 1:n)
    {risk.set[(n-num.risk[i]+1):n,i]=1}
 
 score.mat=t((t(x)<cutp)*1)
 score.tot=apply(score.mat,1,sum)       
 score.bar=rev(cumsum(rev(score.tot)))/(n:1)

 u0.stat=sum(delta*(score.tot-score.bar[time.index]))
 i0.stat=sum(delta/num.risk*((rev(cumsum(rev(score.tot^2)))-score.bar^2*(n:1)))[time.index])
 zscore=u0.stat/sqrt(i0.stat+1e-8)
 

 u.pool=i1.pool=matrix(0, n, p)
 for(i in 1:p)
    {idx=order.index[,i]
     idy=order(idx)

     u.pool[,i]=(cumsum(delta/num.risk))[time.index.plus][idx]

     risk.set.sub=risk.set[idx,]
     w=apply(risk.set.sub,2,cumsum)
     w=rbind(0, w[-n,])
     w=t(w[idy,])

     i1.pool[,i]=(cumsum(delta*(num.risk-1)/num.risk^2))[time.index.plus][idx]-2*(diag((apply(delta/num.risk^2*w,2,cumsum))[time.index.plus,]))[idx]
     }


 count=1
 loop=1
 while(loop==1)
      {score.mat=t((t(x)<cutp)*1)
       score.tot=apply(score.mat,1,sum)       
   
       score.pool=score.bar.pool=matrix(0, n, p)
       for(i in 1:p)
          {idx=order.index[,i] 
           score.pool[,i]=(score.tot-score.mat[,i])
           score.bar.pool[,i]=rev(cumsum(rev(score.tot-score.mat[,i])))/(n:1)
           }
       
       i.pool=matrix(0, n, p)
       for(i in 1:p)
           {
            idx=order.index[,i]
            i.pool[,i]=2*u.pool[,i]*score.pool[idx,i]-2*(cumsum(delta/num.risk*score.bar.pool[time.index,i]))[time.index.plus][idx]
            } 
       i.pool=i.pool+i1.pool

       u0.stat=apply(delta*(score.pool-score.bar.pool[time.index,]),2,sum)
       i0.stat=apply(delta*(1/(n:1)*((apply((score.pool^2)[n:1,],2,cumsum))[n:1,]-score.bar.pool^2*(n:1)))[time.index,],2,sum)
       
    
       score.stat=t(u0.stat+t(apply(delta.pool-u.pool,2,cumsum)))
       v.stat=t(i0.stat+t(apply(i.pool,2,cumsum)))

       for(i in 1:p)
          {idx=x.replicate[,i]
           v.stat[,i]=v.stat[idx ,i]   
           score.stat[,i]=score.stat[idx, i]
          } 

       test=score.stat/sqrt(v.stat+1e-8)

       if(mincut==0)
         {mtest=apply(test,2,max)}

       if(mincut>0)
         {mtest=rep(0, p)
          for(i in 1:p)
             {
              test.post=test[max(effect.range)+1,i]
              test0=test[effect.range,i]
              test0=test0[test0!=test.post]
              mtest[i]=max(test0)
             }
          }
       loop=0
       if(max(mtest)>zscore[count]+1e-8)
         {loop=1
          i=(1:p)[mtest==max(mtest)][1]
          i0=max(effect.range[test[effect.range,i]==mtest[i]])+1
          if(i0>n)
            mcut=Inf
          if(i0<=n)
            mcut=x[order.index[,i],i][i0]
          cutp[i]=mcut
          zscore=c(zscore, max(mtest)) 
          count=count+1           
         }        
       }
       
       
  return(list(cutp=cutp, zscore=zscore))
 }



############################################################################################################################################
############################################################################################################################################
############################################################################################################################################

backfit.cox.interaction=function(x, trt, y, delta, cutp, mincut=0)
{n=length(y)
 p=length(x[1,])
 id=order(y)
 y=y[id]
 delta=delta[id]
 trt=trt[id]
 x=x[id,]

 if(mincut>0)
   {effect.range=ceiling(n*mincut):floor(n*(1-mincut))}
 if(mincut==0)
   {effect.range=1:n}

 ntime=length(unique(y))
 time.index=rep(c(1, cumsum(table(y))+1)[-(ntime+1)], table(y))
 num.risk=(n:1)[time.index]
 time.index.plus=rep(cumsum(table(y)), table(y))


 fit=coxph(Surv(y, delta)~trt, method="breslow")
 beta=fit$coef
 sigma0=fit$var
 risk=exp(beta*trt)

 S0=(rev(cumsum(rev(risk))))[time.index]
 Sr=(rev(cumsum(rev(risk*trt))))[time.index]


 order.index=risk.pool=delta.pool=trt.pool=x.replicate=matrix(0, n, p)
 for(i in 1:p)
    {idx=order(x[,i])
     order.index[,i]=idx
     risk.pool[,i]=risk[idx]
     delta.pool[,i]=delta[idx]
     trt.pool[,i]=trt[idx]
     x.replicate[,i]=rep((1:n)[cumsum(table(x[,i]))], table(x[,i]))
     }

 ds0=(cumsum(delta/S0))[time.index.plus]
 drs0=(cumsum(delta*Sr/S0^2))[time.index.plus]
 dss0=(cumsum(delta/S0^2))[time.index.plus]

 ds0.pool=drs0.pool=dss0.pool=matrix(0, n, p)
 for(i in 1:p)
    {idx=order(x[,i])
     ds0.pool[,i]=ds0[idx]
     dss0.pool[,i]=dss0[idx]
     drs0.pool[,i]=drs0[idx]
     }

 score.mat=t((t(x)<cutp)*1)
 score.tot=apply(score.mat,1,sum)       

 Swwrr0=(rev(cumsum(rev(risk*trt^2*score.tot^2))))[time.index]
 Swrr0=(rev(cumsum(rev(risk*trt^2*score.tot))))[time.index]
 Swr0=(rev(cumsum(rev(risk*trt*score.tot))))[time.index]

 i1=sum(delta*Swwrr0/S0)
 i2=sum(delta*Swr0^2/S0^2)
 i3=sum(delta*Swrr0/S0)
 i4=sum(delta*Swr0*Sr/S0^2)
 v.stat=(i1-i2)-(i3-i4)^2*sigma0[1,1]
 score.stat=sum(delta*score.tot*trt)-sum(delta*Swr0/S0)
 zscore=score.stat/sqrt(v.stat+1e-8)
 
 
 subject=matrix(0, n, n)
 for(i in 1:n)
    {subject[(n-num.risk[i]+1):n,i]=1}
 
 i2.add=apply(trt.pool^2*risk.pool^2*dss0.pool,2,cumsum)
 i3.pool=apply(trt.pool^2*risk.pool*ds0.pool,2,cumsum)
 i4.pool=apply(trt.pool*risk.pool*drs0.pool,2,cumsum)

 w1.pool=matrix(0, n, p)
 for(i in 1:p)
    {idx=order.index[,i]
     subject.id=subject[idx,]
     temp1=apply(risk[idx]*trt[idx]*subject.id,2,cumsum)
     temp1=rbind(0, temp1[-n,])    
     w1.pool[,i]=diag((apply(t(temp1)*delta/S0^2,2,cumsum))[time.index.plus,][idx,])
     }

 score.stat.pool=apply(delta.pool*trt.pool,2,cumsum)-apply(trt.pool*risk.pool*ds0.pool,2,cumsum)

 count=1
 loop=1
 while(loop==1)
       {score.mat=t((t(x)<cutp)*1)
        score.tot=apply(score.mat,1,sum)       
   
             
        score.pool=score.pool0=matrix(0, n, p)
        for(i in 1:p)
           {idx=order.index[,i]
            score.pool[,i]=(score.tot-score.mat[,i])[idx] 
            score.pool0[,i]=(score.tot-score.mat[,i])}


        Swwrr0=((apply((risk*trt^2*score.pool0^2)[n:1,],2,cumsum))[n:1,])[time.index,]
        Swrr0=((apply((risk*trt^2*score.pool0)[n:1,],2,cumsum))[n:1,])[time.index,]
        Swr0=((apply((risk*trt*score.pool0)[n:1,],2,cumsum))[n:1,])[time.index,]


        w.pool=matrix(0, n, p)
        for(i in 1:p)
           {idx=order.index[,i]  
            temp2=(cumsum((rev(cumsum(rev(risk*trt*score.pool0[,i]))))[time.index]*delta/S0^2))[time.index.plus]
            w.pool[,i]=w1.pool[,i]+temp2[idx]
            }

        i10=apply(delta*Swwrr0/S0,2,sum)
        i20=apply(delta*Swr0^2/S0^2,2,sum)
        i30=apply(delta*Swrr0/S0,2,sum)
        i40=apply(delta*Swr0*Sr/S0^2,2,sum)   

        i1.pool=apply((2*score.pool+1)*trt.pool^2*risk.pool*ds0.pool,2,cumsum)
        i2.pool=apply(2*trt.pool*risk.pool*w.pool,2,cumsum)+i2.add
       
        i1=t(t(i1.pool)+i10)
        i2=t(t(i2.pool)+i20)
        i3=t(t(i3.pool)+i30)
        i4=t(t(i4.pool)+i40)

        v.stat=(i1-i2)-(i3-i4)^2*sigma0[1,1]

        score.stat0=apply(delta.pool*score.pool*trt.pool,2,sum)-apply(delta*Swr0/S0,2,sum)
        score.stat=t(t(score.stat.pool)+score.stat0)

        for(i in 1:p)
            {idx=x.replicate[,i]
             v.stat[,i]=v.stat[idx ,i]   
             score.stat[,i]=score.stat[idx, i]
            } 

         test=score.stat/sqrt(v.stat+1e-8)

         if(mincut==0)
           {mtest=apply(test,2,max)}

         if(mincut>0)
           {mtest=rep(0, p)
            for(i in 1:p)
               {
                test.post=test[max(effect.range)+1,i]
                test0=test[effect.range,i]
                test0=test0[test0!=test.post]
                mtest[i]=max(test0)
             }
            }
         loop=0
         if(max(mtest)>zscore[count]+1e-8)
           {loop=1
            i=(1:p)[mtest==max(mtest)][1]
            i0=max(effect.range[test[effect.range,i]==mtest[i]])+1
            if(i0>n)
               mcut=Inf
            if(i0<=n)
              mcut=x[order.index[,i],i][i0]
            cutp[i]=mcut
            zscore=c(zscore, max(mtest)) 
            count=count+1           
            }        
          }
       
       
     return(list(cutp=cutp, zscore=zscore))
    }
