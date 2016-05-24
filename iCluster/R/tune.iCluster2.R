tune.iCluster2=function(datasets, k, n.lambda=NULL, nrep=10, mc.cores=6, max.iter=10)
  {
    n <-dim(datasets[[1]])[1]
    m=length(datasets) 
    p=unlist(lapply(1:m,function(l){dim(datasets[[l]])[2]}))
    sum.p=sum(p)

    #pre-checks
    if(is.list(datasets)==F)stop("datasets must be a list")
    if(sum(is.na(datasets))>0)stop("data cannot have NAs. please exclude or impute missing values.")
    if(missing(k))stop("must specify the number of clusters k")
    
    #uniform sampling of penalty parameter values (Fang and Wang 1995)
    data(glp)
    n.glp=c(10,21,35,307,1069)
    lbd = list()
    lbd[[1]] = 1:1000
    lbd[[2]] = c(8, 13, 21, 34, 55, 89, 144, 233, 377, 610)
    lbd[[3]] = c(35, 101, 135, 185, 266, 418, 579, 828, 1010)
    lbd[[4]] = c(307, 526, 701, 1019, 2129, 3001, 4001, 5003, 6007)
    lbd[[5]] = c(1069,1543, 2129, 3001, 4001, 5003, 6007, 8191)
    
    if(is.null(n.lambda)){
      n.lambda=n.glp[m]
      cat(n.glp[m]," points of lambdas are sampled by uniform design.\n")
    }else{
      if(any(n.lambda==lbd[[m]])){
        cat(n.lambda," points of lambdas are sampled by uniform design.\n")
      }else if(m>1){
        cat("Error: n.lambda is not a valid value\n")
        cat("Valid value for n.lambda are ",lbd[[m]],"\n")
        stop()
      }else{
        warning()
        cat(n.glp[m]," points of lambdas are sampled by uniform design; out of range.\n")
      }
    }
    
    which.glp=paste("s",m, ".txt", sep='')
    h=glp[[which.glp]]
    h=h[-1,which(h[1,]==n.lambda)]
    h=c(1,h)
    if(m==1){ud=as.matrix(seq(1,2*n.lambda-1, by=2)/2/n.lambda,nrow=n.lambda)}
    if(m>1){
      ud=matrix(NA,nrow=n.lambda, ncol=m)
      for(s in 1:m){
        ud[,s]=((1:n.lambda)*h[s]-0.5)/n.lambda  
      }
    }
    ud=ud-floor(ud)
        
    
    
    cat("Begin parallel computation... \n")
    fit = mclapply(1:nrow(ud), 
                   FUN=function(x)iCluster.RI(datasets,lambda=ud[x,],k,nrep,
                         max.iter), mc.silent=T, mc.cores=mc.cores, mc.preschedule =F)    
    RI=unlist(fit)
    best.fit=iCluster2(datasets=datasets,k=k, lambda=ud[which.max(RI),], verbose=F, max.iter=max.iter)   
    out=list(best.fit=best.fit, ud=ud, RI=RI)
    return(out)
   
}    

iCluster.RI=function(datasets, k, lambda, max.iter,nrep){

    n <-dim(datasets[[1]])[1]
    m=length(datasets) 
    p=unlist(lapply(1:m,function(l){dim(datasets[[l]])[2]}))
    sum.p=sum(p)
      
    ps.c=ps.adjusted.c=NULL    
    for(nsplit in 1:nrep){
      
      folds <- split(sample(seq(n)), rep(1:2, length=n))
      
      omit <- folds[[1]]
      trdata=alist()
      tsdata=alist()
      for(j in 1:m){
        trdata[[j]] <- datasets[[j]][-omit, ]
        tsdata[[j]] <- datasets[[j]][omit, ]
      }
      fit.trdata=iCluster2(datasets=trdata, k=k, lambda=lambda,verbose=F, max.iter=max.iter)
      W=fit.trdata$W
      PSI=fit.trdata$PSI
      sigma=W%*%t(W)+diag(PSI)
      stacked.tsdata=matrix(NA, nrow=length(omit), ncol=sum.p)
      for (t in 1:m) {
        if(t==1){idx=1:p[t]}else{idx=(sum(p[1:(t-1)])+1):sum(p[1:t])}
        stacked.tsdata[,idx] <- tsdata[[t]]
      }
      pred.expZ=t(W)%*%solve(sigma)%*%t(stacked.tsdata)
      
      
      fit=iCluster2(datasets=tsdata,k=k, lambda=lambda,verbose=F, max.iter=max.iter)    
      pred=predict.kmeans(fit.trdata,t(pred.expZ))
      ps.adjusted.c[nsplit]=aRandIndex(fit$clusters, pred) 
      
    }
    
    ps.adjusted=median(ps.adjusted.c)
    return(ps.adjusted)
  
}

aRandIndex=function (x, y) 
{
    x <- as.vector(x)
    y <- as.vector(y)
    xx <- outer(x, x, "==")
    yy <- outer(y, y, "==")
    upper <- row(xx) < col(xx)
    xx <- xx[upper]
    yy <- yy[upper]
    a <- sum(as.numeric(xx & yy))
    b <- sum(as.numeric(xx & !yy))
    c <- sum(as.numeric(!xx & yy))
    d <- sum(as.numeric(!xx & !yy))
    ni <- (b + a)
    nj <- (c + a)
    abcd <- a + b + c + d
    q <- (ni * nj)/abcd
    (a - q)/((ni + nj)/2 - q)
}
RandIndex=function (c1, c2)
{
    c1 <- as.vector(c1)
    c2 <- as.vector(c2)
    xx <- outer(c1, c1, "==")
    yy <- outer(c2, c2, "==")
    upper <- row(xx) < col(xx)
    xx <- xx[upper]
    yy <- yy[upper]
    a <- sum(as.numeric(xx & yy))
    d <- sum(as.numeric(!xx & !yy))
    (a+d)/choose(length(c2),2)
}

predict.kmeans=function(km, data)
{k <- nrow(km$centers)
n <- nrow(data)
d <- as.matrix(dist(rbind(km$centers, data)))[-(1:k),1:k]
out <- apply(d, 1, which.min)
return(out)}
