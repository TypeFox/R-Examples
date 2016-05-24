cv.lm.main=function(x, y, K.cv=5, num.replicate=1, nsteps, mincut=0.1, backfit=F, maxnumcut=1, dirp=0){
x=as.matrix(x)
n=length(y)

sc.tot=matrix(0, nrow=K.cv, ncol=nsteps)
preval.tot=matrix(0, nrow=nsteps, ncol=n)

for(b in 1:num.replicate)
   {g=sample(rep(1:K.cv,(n+K.cv-1)/K.cv))
    g=g[1:n]
    cva=vector("list",K.cv)

    sc=matrix(0, nrow=K.cv, ncol=nsteps)
    preval=matrix(0, nrow=nsteps, ncol=n)
    for(i in 1:K.cv)
       {cva[[i]]=lm.main(x[g!=i ,],  y[g!=i], mincut=mincut, nsteps=nsteps, backfit=backfit, maxnumcut=maxnumcut, dirp=dirp)
        for(ii in 1:nsteps)
           {aa=index.prediction(cva[[i]]$res[[ii]],x[g==i,])
            fit=lm(y[g==i]~aa)
            sc[i,ii]=summary(fit)$coef[2,3]
            preval[ii,g==i]=aa
            }
        }
    sc.tot=abs(sc.tot)+sc
    preval.tot=preval.tot+preval
    }

meansc=colMeans(sc.tot)
kmax=which.max(meansc)
pvfit.score=rep(0, nsteps)
for(i in 1:nsteps)
   {fit=lm(y~preval.tot[i,])
    pvfit.score[i]=summary(fit)$coef[2,3]
    }


return(list(kmax=kmax, meanscore=meansc/num.replicate, pvfit.score=pvfit.score, preval=preval.tot/num.replicate))
}

