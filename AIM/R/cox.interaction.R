cox.interaction=function(x, trt, y, delta, nsteps=8, mincut=0.1, backfit=F, maxnumcut=1, dirp=0)
{
n=length(y)
id=order(y)
y=y[id]
x=x[id,]
delta=delta[id]
trt=trt[id]

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

ntime=length(unique(y))
time.index=rep(c(1, cumsum(table(y))+1)[-(ntime+1)], table(y))
num.risk=(n:1)[time.index]
time.index.plus=rep(cumsum(table(y)), table(y))


id.in=NULL
id.out=1:p.true
p.out=p
zvalue=NULL
cut.value=NULL
direction=NULL
imax=NULL
res=as.list(NULL)
num.cut=rep(0, p.true)

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


subject=matrix(0, n, n)
for(i in 1:n)
    {subject[(n-num.risk[i]+1):n,i]=1}


i2.diff.pool=apply(trt.pool^2*risk.pool^2*dss0.pool,2,cumsum)
i3.diff.pool=apply(trt.pool^2*risk.pool*ds0.pool,2,cumsum)
i4.diff.pool=apply(trt.pool*risk.pool*drs0.pool,2,cumsum)
score.diff.pool=apply(delta.pool*trt.pool,2,cumsum)-apply(trt.pool*risk.pool*ds0.pool,2,cumsum)

w1.pool=matrix(0, n, p.out)
for(i in 1:p.out)
   {idx=order.index[,i]
    subject.id=subject[idx,]
    temp1=apply(risk[idx]*trt[idx]*subject.id,2,cumsum)
    temp1=rbind(0, temp1[-n,])
    w1.pool[,i]=diag((apply(t(temp1)*delta/S0^2,2,cumsum))[time.index.plus,][idx,])
    }
     
count=1
loop=1
while(loop==1){


Swwrr0=(rev(cumsum(rev(risk*trt^2*score.current^2))))[time.index]
Swrr0=(rev(cumsum(rev(risk*trt^2*score.current))))[time.index]
Swr0=(rev(cumsum(rev(risk*trt*score.current))))[time.index]


score.pool=w.pool=matrix(0, n, p.out)
for(i in 1:p.out)
    {idx=order.index[,i]   
     score.pool[,i]=score.current[idx]
     temp=(cumsum((rev(cumsum(rev(risk*trt*score.current))))[time.index]*delta/S0^2))[time.index.plus]
     w.pool[,i]=w1.pool[,i]+temp[idx]
     }
     

i1=sum(delta*Swwrr0/S0)+apply((2*score.pool+1)*trt.pool^2*risk.pool*ds0.pool,2,cumsum)
i2=sum(delta*Swr0^2/S0^2)+apply(2*trt.pool*risk.pool*w.pool,2,cumsum)+i2.diff.pool
i3=sum(delta*Swrr0/S0)+i3.diff.pool
i4=sum(delta*Swr0*Sr/S0^2)+i4.diff.pool


v.stat=(i1-i2)-(i3-i4)^2*sigma0[1,1]
score.stat=sum(delta*score.current*trt)-sum(delta*Swr0/S0)+score.diff.pool

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
         

i.sel=which.max(mtest)
score.current=score.current+(x.out[,i.sel]<mcut[i.sel])

marker.sel=marker.x[i.sel]
id.in=c(id.in, marker.sel)
direction=c(direction, direction.x[i.sel])
cut.value=c(cut.value, -mcut[i.sel]*direction.x[i.sel])
num.cut[marker.sel]=num.cut[marker.sel]+1


if(num.cut[marker.sel]==maxnumcut)
  {id.exclude=(1:p.out)[marker.x==marker.sel]

   x.out=x.out[,-id.exclude, drop=F]
   order.index=order.index[, -id.exclude, drop=F]
   trt.pool=trt.pool[, -id.exclude, drop=F]
   risk.pool=risk.pool[, -id.exclude, drop=F]
   ds0.pool=ds0.pool[, -id.exclude, drop=F]
   i2.diff.pool=i2.diff.pool[, -id.exclude, drop=F]
   i3.diff.pool=i3.diff.pool[, -id.exclude, drop=F]
   i4.diff.pool=i4.diff.pool[, -id.exclude, drop=F]
   score.diff.pool=score.diff.pool[, -id.exclude, drop=F]
   w1.pool=w1.pool[, -id.exclude, drop=F]
   x.replicate=x.replicate[, -id.exclude, drop=F]
   
   direction.x=direction.x[-id.exclude]
   marker.x=marker.x[-id.exclude]   
   p.out=p.out-length(id.exclude)                  
   }

       
if(backfit==T)
  {if(count>1)
     {x.adj=x0[, id.in]
      x.adj=-t(t(x.adj)*direction)
      cutp=-cut.value*direction
      fit=backfit.cox.interaction(x.adj, trt, y, delta, cutp, mincut=mincut)
         
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
