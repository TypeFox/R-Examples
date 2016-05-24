
logistic.interaction=function(x, trt, y, nsteps=8, mincut=0.1, backfit=F, maxnumcut=1, dirp=0, weight=1)
# interaction AIM for logistic
#  y- n vector of outcomes (0-1)
# trt= 0-1 treatment variable- should be centered before input into function
# x- n by p matrix of predcitors
# nsteps= # number terms to add 
# mincut- min  cutting prop (at either end)  
# backfit- should split points be re-adjusted after every term is added?
# maxnumcut- max # of splits per feature
# dirp- directions of splits allowed- a p-vector of 0,1, or -1.
#  1 means split to the left; -1 to the right; 0 means either. dirp=0
#   means all splits are unrestricted
#  weight= obs wt given to second class

# Output is a list res, with each component summarizing the model
# after a term has been added. columns are
#  jmax- predictor#,  cutp- cutpoint,  maxdir- direction of split (1 means le,
# -1 means >) and    maxsc- likelihood score achieved
{n=length(y)
 x0=x
 p.true=length(x0[1,])

 if(length(weight)==1)
   {wt=y*(weight-1)+1}
 if(length(weight)==n)
   {wt=weight} 

 x=cbind(x0, -x0)
 marker.x=rep(1:p.true, 2)
 direction.x=c(rep(-1, p.true), rep(1, p.true))

 if(sum(abs(dirp))>0)
    {index1=(1:p.true)[dirp==-1 | dirp==0]
     index2=(1:p.true)[dirp==1 | dirp==0]
     x=cbind(x0[,index1], -x0[,index2])
     marker.x=c((1:p.true)[index1], (1:p.true)[index2])
     direction.x=c(rep(-1, length(index1)), rep(1, length(index2)))
     }                   
 p=length(x[1,])
 

 if(mincut>0)
   {effect.range=ceiling(n*mincut):floor(n*(1-mincut))}
 if(mincut==0)
   {effect.range=1:n}

 score.current=rep(0,n)
 num.cut=rep(0, p.true)

 order.index=matrix(0, n, p)
 for(i in 1:p.true)
    {order.index[,i]=order(x[,i])
     order.index[,i+p.true]=rev(order.index[,i])}

 wt.pool=trt.pool=y.pool=x.replicate=matrix(0, n, p)
 for(i in 1:p)
    {y.pool[,i]=y[order.index[,i]]
     trt.pool[,i]=trt[order.index[,i]]
     wt.pool[,i]=wt[order.index[,i]]
     x.replicate[,i]=rep((1:n)[cumsum(table(x[,i]))], table(x[,i]))}

 x.out=x 
 id.in=NULL
 p.out=p
 zvalue=NULL
 cut.value=NULL
 direction=NULL
 res=as.list(NULL)

 fit=glm(y~trt, family="binomial", weight=wt)
 phat=fit$fitt
 sigma=sum(summary(fit)$cov.scale)


 count=1
 loop=1
 while(loop==1)
      {score.pool=phat.pool=matrix(0, n, p.out)
       for(i in 1:p.out)
          {score.pool[,i]=score.current[order.index[,i]]
           phat.pool[,i]=phat[order.index[,i]]}

       score.stat=apply(trt.pool*wt.pool*(y.pool-phat.pool),2,cumsum)+sum((y-phat)*trt*score.current*wt)
       v11.stat=2*apply(score.pool*trt.pool*wt.pool^2*phat.pool*(1-phat.pool),2,cumsum)+apply(trt.pool*wt.pool^2*phat.pool*(1-phat.pool),2,cumsum)+sum((score.current*wt)^2*trt*phat*(1-phat))
       v22.stat=(sum(score.current*trt*wt*phat*(1-phat))+apply(trt.pool*wt.pool*phat.pool*(1-phat.pool),2,cumsum))^2*sigma
       v.stat=(v11.stat-v22.stat)

       
       for(i in 1:p.out)
          {idx=x.replicate[,i]
           v.stat[,i]=v.stat[idx ,i]   
           score.stat[,i]=score.stat[idx, i]
          } 

       test=score.stat/sqrt(v.stat+1e-8)

       if(mincut==0)
         {mtest=apply(test,2,max)}

       if(mincut>0)
         {mtest=rep(0, p.out)
          for(i in 1:p.out)
             {
              test.post=test[max(effect.range)+1,i]
              test0=test[effect.range,i]
              test0=test0[test0!=test.post]
              mtest[i]=max(test0)}
          }


       mcut=rep(0, p.out)
       i0=rep(0, p.out)
       for(i in 1:p.out)
          {i0[i]=max(effect.range[test[effect.range,i]==mtest[i]])+1
           if(i0[i]>n)
              mcut[i]=Inf
           if(i0[i]<=n)
              mcut[i]=x.out[order.index[,i],i][i0[i]]
           }
         
       i.sel=(1:p.out)[mtest==max(mtest)][1]
       score.current=score.current+(x.out[,i.sel]<mcut[i.sel])


       marker.sel=marker.x[i.sel]
       id.in=c(id.in, marker.sel)
       direction=c(direction, direction.x[i.sel])
       cut.value=c(cut.value, -mcut[i.sel]*direction.x[i.sel])
       num.cut[marker.sel]=num.cut[marker.sel]+1

       if(num.cut[marker.sel]==maxnumcut)
         {id.exclude=(1:p.out)[marker.x==marker.sel]

          x.out=x.out[,-id.exclude,drop=F]
          order.index=order.index[, -id.exclude,drop=F]
          y.pool=y.pool[,-id.exclude,drop=F]
          wt.pool=wt.pool[,-id.exclude, drop=F]
          trt.pool=trt.pool[,-id.exclude,drop=F]
          x.replicate=x.replicate[, -id.exclude,drop=F]          

          direction.x=direction.x[-id.exclude]
          marker.x=marker.x[-id.exclude]
          p.out=p.out-length(id.exclude)
         }

       
       if(backfit==T)
         {if(count>1)
            {x.adj=x0[, id.in]
             x.adj=-t(t(x.adj)*direction)
             cutp=-cut.value*direction
             fit=backfit.logistic.interaction(x.adj, trt, y, cutp, mincut=mincut, wt=wt)
         
             jmax=id.in
             cutp=-fit$cutp*direction
             maxdir=direction
             
             res[[count]]=cbind(jmax, cutp, maxdir)
             zvalue=c(zvalue,  max(fit$zscore)) 
             score.current=apply((t(x.adj)<fit$cutp),2,sum)
             }

          if(count==1)
            {jmax=id.in
             cutp=cut.value
             maxdir=direction
             res[[count]]=cbind(jmax, cutp, maxdir)
             zvalue=c(zvalue,  max(mtest))
             }
           }
       if(backfit==F)
         {cutp=cut.value
          maxdir=direction
          jmax=id.in
          zvalue=c(zvalue, max(mtest))
          maxsc=zvalue
          res[[count]]=cbind(jmax, cutp, maxdir, maxsc)
          }
       
       count=count+1
       loop=(length(id.in)<nsteps)
       }

       
 return(list(res=res, maxsc=zvalue))
 }

