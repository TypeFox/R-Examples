
cox.main=function(x, y, delta, nsteps=8, mincut=0.1, backfit=F, maxnumcut=1, dirp=0)
{
 n=length(x[,1])
 id=order(y)
 y=y[id]
 delta=delta[id]
 x=x[id,]

 
 x0=x
 p.true=length(x0[1,])
 
 x=cbind(x0, -x0)
 marker.x=rep(1:p.true, 2)
 direction.x=c(rep(-1, p.true), rep(1, p.true))

 if(sum(abs(dirp))>0)
    {index1=(1:p.true)[dirp==-1 | dirp==0]
     index2=(1:p.true)[dirp==1 | dirp==0]
     x=cbind(x0[,index1], x0[,index2])
     marker.x=c((1:p.true)[index1], (1:p.true)[index2])
     direction.x=c(rep(-1, length(index1)), rep(1, length(index2)))
     }                   
 p=length(x[1,])
 

 if(mincut>0)
   {effect.range=ceiling(n*mincut):floor(n*(1-mincut))}
 if(mincut==0)
   {effect.range=1:n}

 score.current=rep(0,n)
 x.out=x

 
 order.index=delta.pool=x.replicate=matrix(0, n, p)
 for(i in 1:p)
      {idx=order(x[,i])
       order.index[,i]=idx
       delta.pool[,i]=delta[idx]
       x.replicate[,i]=rep((1:n)[cumsum(table(x[,i]))], table(x[,i]))
       } 
 
 id.in=NULL
 id.out=1:p.true
 p.out=p
 zvalue=NULL
 cut.value=NULL
 direction=NULL
 imax=NULL
 res=as.list(NULL)
 num.cut=rep(0, p.true)

 tby=table(y)
 ntime=length(unique(y))
 time.index=rep(c(1, cumsum(tby)+1)[-(ntime+1)], tby)
 num.risk=(n:1)[time.index]
 time.index.plus=rep(cumsum(tby), tby)


 risk.set=matrix(0, n, n)
 for(i in 1:n)
    {risk.set[(n-num.risk[i]+1):n,i]=1}


 
 i1.pool=i3.pool=u.pool=matrix(0, n, p.out)
 u0=(cumsum(delta/num.risk))[time.index.plus]
 i10=(cumsum(delta*(num.risk-1)/num.risk^2))[time.index.plus]
 for(i in 1:p.out)
    {idx=order.index[,i]
     idy=order(idx)

     u.pool[,i]=u0[idx]

     i1.pool[,i]=i10[idx]
      
     risk.set.sub=risk.set[idx,]
     w=apply(risk.set.sub,2,cumsum)
     w=rbind(0, w[-n,])
     w=t(w[idy,])
     i3.pool[,i]=2*(diag((apply(delta/num.risk^2*w,2,cumsum))[time.index.plus,]))[idx]
     }

 
 count=1
 loop=1
 while(loop==1)
  {score.bar=(rev(cumsum(rev(score.current)))/(n:1))
   u0.stat=sum(delta*(score.current-score.bar[time.index]))
   i0.stat=sum(delta/num.risk*((rev(cumsum(rev(score.current^2)))-score.bar^2*(n:1)))[time.index])

 
   i.pool=matrix(0, n, p.out)
   for(i in 1:p.out)
      {
       idx=order.index[,i]                
       i.pool[,i]=i1.pool[,i]+2*u.pool[,i]*score.current[idx]-2*(cumsum(delta/num.risk*score.bar[time.index]))[time.index.plus][idx]-i3.pool[,i]
       } 


     score.stat=u0.stat+apply(delta.pool-u.pool,2,cumsum)
     v.stat=i0.stat+apply(i.pool,2,cumsum)

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
             {test.post=test[max(effect.range)+1,i]
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
        delta.pool=delta.pool[,-id.exclude,drop=F]
        x.replicate=x.replicate[, -id.exclude,drop=F]
        i1.pool=i1.pool[,-id.exclude,drop=F]
        i3.pool=i3.pool[,-id.exclude,drop=F]
        u.pool=u.pool[,-id.exclude,drop=F]
        direction.x=direction.x[-id.exclude]
        marker.x=marker.x[-id.exclude]   
        p.out=p.out-length(id.exclude)                  
        }
          

       
     if(backfit==T)
         {if(count>1)
            {x.adj=x0[, id.in]
             x.adj=-t(t(x.adj)*direction)
             cutp=-cut.value*direction
             fit=backfit.cox.main(x.adj, y, delta, cutp, mincut=mincut)
         
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
