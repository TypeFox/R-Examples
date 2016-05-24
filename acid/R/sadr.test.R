sadr.test <-
  function(data,y.pos=NULL,dist1,dist2,params.m,mcmc=TRUE,mcmc.params.a,ygrid, bsrep=10,n.startvals=300,dist.para.table) {
    helpfun1 <- function(cdfobj,y,ygrid,...) {approx(x=ygrid,y=sort(cdfobj),xout=y,rule=2,method="linear")$y} #x pos. on grid, y density at x, xout-where to give interpolated density
    
    if(is.null(y.pos)) y.pos<-1
    n <- dim(data)[1] 
    if(n!=dim(params.m)[1]) print("Data not of same length as parameter matrix!")
    pn1 <- dist.para.table[which(dist.para.table$dist == dist1),3] #number of parameters dist1
    pn2 <- dist.para.table[which(dist.para.table$dist == dist2),3] #number of parameters dist2
    pn  <- pn1+pn2 #number of parameters both dist
    
    cond.cdf<-matrix(NA,length(ygrid),n)
    for(i in 1:n){
      cond.cdf[,i]<-pval.md(ygrid,dist1,dist2,
                            theta=params.m[i,1:pn],
                            params.m[i,pn+1],params.m[i,pn+2],params.m[i,pn+3],dist.para.table)
    }
    
    y<-as.matrix(data[,y.pos])
    x<-as.matrix(data[,-c(y.pos)])
    
    startvals <- sort(sample(1:n,round(min(c(n, n.startvals))),replace=FALSE)) #maximal 300 Startwerte to estimate the test statistic - sufficient??
    x.ev <- as.matrix(x[startvals,])
    y.ev <- y[startvals] 
    
    cond.cdf.X <- apply(cond.cdf,2, helpfun1,y=y.ev,ygrid=ygrid)  #row i gives cdf values for individual i, as interpolated using the cdf given X values of individual j, d.h. Fhat_n(y|X-i)
    
    cdfdiff <- rep(NA,length(y.ev))
    for (i in 1:length(y.ev)) {
      cdfdiff[i] <- mean(( cond.cdf.X[i,]  - I(y <= y.ev[i]))*apply(t(x)- x.ev[i,] <= 0, 2,prod)) #1/n*(FN(y|X_i)-I(Y_i <= y))*I(X_i <= x)
    }
    
    teststat.ks  <- sqrt(n)*max(abs(cdfdiff)) #ks measure
    teststat.cvm <- n*sum(cdfdiff^2) # equivalent to the sum, i.e. T_n
    
    ## Start the bootstrap procedure
    teststat.ks.bs  <- rep(NA,bsrep)
    teststat.cvm.bs <- rep(NA,bsrep)
    
    for (b in 1:bsrep) {
      if(mcmc==TRUE){
        m<-sample(1:dim(mcmc.params.a)[1],1) #sampling MCMC draw for all pi parameters
        p0.bs<-mcmc.params.a[m,,pn+1] #sampling MCMC draw for p0
        p1.bs<-mcmc.params.a[m,,pn+4] #sampling MCMC draw for p1
        m<-sample(1:dim(mcmc.params.a)[1],1) #sampling MCMC draw for all ft parameters
        theta1.bs<-mcmc.params.a[m,,1:pn1]
        m<-sample(1:dim(mcmc.params.a)[1],1) #sampling MCMC draw for all pt parameters
        theta2.bs<-mcmc.params.a[m,,(pn1+1):pn]
        theta.bs<-cbind(theta1.bs,theta2.bs)
      }else{
        p0.bs<-params.m[,pn+1] #sampling MCMC draw for p0
        p1.bs<-params.m[,pn+4] #sampling MCMC draw for p1
        theta.bs<-params.m[,1:pn]
      }
      
      sam.bs<-sample(1:n,n,replace=TRUE) 
      data.bs <- data[sam.bs,] 
      x.bs   <- as.matrix(x[sam.bs,])
      y.bs   <- matrix(NA,dim(x.bs)[1],1) 
      for(i in 1:n){
        si<-sam.bs[i]
        p0<-p0.bs[si]
        ppt<-(1-p0.bs[si])*(1-p1.bs[si])
        pft<-(1-p0.bs[si])*(p1.bs[si])
        data.bs[i,y.pos]<-y.bs[i,1]<-ysample.md(1,dist1,dist2,
                                                theta=theta.bs[si,],p0,pft,ppt,dist.para.table)
      }
      
      startvals <- sort(sample(1:n,round(min(c(n, n.startvals))),replace=FALSE))
      x.ev.bs <- as.matrix(x.bs[startvals,])
      y.ev.bs <- y.bs[startvals] 
      
      cond.cdf.bs<-matrix(NA,length(ygrid),n)
      for(i in 1:n){
        si<-sam.bs[i]
        p0<-p0.bs[si] 
        ppt<-(1-p0.bs[si])*(1-p1.bs[si])
        pft<-(1-p0.bs[si])*(p1.bs[si])
        cond.cdf.bs[,i]<-pval.md(ygrid,dist1,dist2,
                                 theta=theta.bs[si,],p0,pft,ppt,dist.para.table)
      }
      cond.cdf.X.bs <- apply(cond.cdf.bs,2, helpfun1,y=y.ev.bs,ygrid=ygrid) 
      
      cdfdiff.bs <- rep(NA,length(y.ev.bs))
      for (i in 1:length(y.ev.bs)) {
        cdfdiff.bs[i] <- mean(( cond.cdf.X.bs[i,]  - I(y.bs <= y.ev.bs[i]))*apply(t(x.bs)- x.ev.bs[i,] <= 0, 2,prod)) #1/n*(FN(y|X_i)-I(Y_i <= y))*I(X_i <= x)
      }
      
      teststat.ks.bs[b] <- sqrt(n)*max(abs(cdfdiff.bs)) 
      teststat.cvm.bs[b] <- n*sum(cdfdiff.bs^2)
    }
    
    pval.ks <- mean(teststat.ks.bs > teststat.ks, na.rm=TRUE)
    pval.cvm <- mean(teststat.cvm.bs > teststat.cvm, na.rm=TRUE)
    
    list(teststat.ks=teststat.ks, pval.ks=pval.ks,teststat.cvm=teststat.cvm, pval.cvm=pval.cvm, test="Structured Additive Distributional Regression",param.distributions=c(dist1,dist2),
         teststat.ks.bs=sort(teststat.ks.bs),teststat.cvm.bs=sort(teststat.cvm.bs))
    
  }
