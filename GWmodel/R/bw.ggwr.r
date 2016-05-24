##Bandwidth selection for generalized GWR models
#Author: Binbin Lu
bw.ggwr<-function(formula, data, family ="poisson", approach="CV",kernel="bisquare",adaptive=FALSE, p=2, theta=0, longlat=F,dMat)
{
	if (is(data, "Spatial"))
  {
    dp.locat<-coordinates(data)
    data <- as(data, "data.frame")
  }
  else
  {
       stop("Given regression data must be Spatial*DataFrame")
  }
  #cat("This selection has been optimised by golden selection.\n")
  mf<- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0L)

  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  y <- model.extract(mf, "response")
  x <- model.matrix(mt, mf)
  dp.n<-nrow(data)
  ####################################Coffee
  if(dp.n>1500)
  {
    cat("Take a cup of tea and have a break, it will take a few minutes.\n")
    cat("          -----A kind suggestion from GWmodel development group\n")
  }
  #################### Recommond to specify a distance matrix
  if (missing(dMat))
      DM.given<-F
  else
  {
    DM.given<-T
    dim.dMat<-dim(dMat)
    if (dim.dMat[1]!=dp.n||dim.dMat[2]!=dp.n)
    stop ("Dimensions of dMat are not correct")
  }
  #########Find the range of the fixed bandwidth
  if(adaptive)
  {
    upper<-dp.n
    lower<-20
  }
  else
  {
    if(DM.given)
    {
      upper<-range(dMat)[2]
      lower<-upper/5000
    }
    ###!!!!!!!! Important note: if the distance matrix is not specified, the running time will be consuming very much by choosing the range of fixed bandwidth when the p is not 2; 
    ### because the range can't be decided by boundary box
    else
    {
      dMat<-NULL
      if (p==2)
      {
        b.box<-bbox(dp.locat)
        upper<-sqrt((b.box[1,2]-b.box[1,1])^2+(b.box[2,2]-b.box[2,1])^2)
        lower<-upper/5000
      }
      else
      {
        upper<-0
        for (i in 1:dp.n)
        {
          dist.vi<-gw.dist(dp.locat=dp.locat, focus=i, p=p, theta=theta, longlat=longlat)
          upper<-max(upper, range(dist.vi)[2])
        }
        lower<-upper/5000
      }
    }
     
  }
  ########################## Now the problem for the golden selection is too computationally heavy
    #Select the bandwidth by golden selection
    bw<-NA    
    if(approach=="cv"||approach=="CV")
       bw <- gold(ggwr.cv,lower,upper,adapt.bw=adaptive,x,y,family=family,kernel,adaptive, dp.locat, p, theta, longlat,dMat)
    else if(approach=="aic"||approach=="AIC"||approach=="AICc")
       bw<-gold(ggwr.aic,lower,upper,adapt.bw=adaptive,x,y,family=family,kernel,adaptive, dp.locat, p, theta, longlat,dMat)    
   # bw<-NA
#    if(approach=="cv"||approach=="CV")
#       bw <- optimize(bw.cv,lower=lower,upper=upper,maximum=FALSE,X=x,Y=y,kernel=kernel,
#       adaptive=adaptive, dp.locat=dp.locat, p=p, theta=theta, longlat=longlat,dMat=dMat,tol=.Machine$double.eps^0.25)
#    else if(approach=="aic"||approach=="AIC"||approach=="AICc")
#       bw<-optimize(bw.aic,lower=lower,upper=upper,x,y,kernel,adaptive, dp.locat, p, theta, longlat,dMat)    
    bw  
}
####Calculate the CV score with a given bandwidth
##Author: Binbin Lu
ggwr.cv<-function(bw, X, Y,family="poisson", kernel="bisquare",adaptive=F, dp.locat,  p=2, theta=0, longlat=F,dMat)
{
   dp.n<-length(dp.locat[,1])
   #########Distance matrix is given or not

  if (is.null(dMat))
      DM.given<-F
  else
  {
    DM.given<-T
    dim.dMat<-dim(dMat)
    if (dim.dMat[1]!=dp.n||dim.dMat[2]!=dp.n)
    stop ("Dimensions of dMat are not correct")
  }
  ############################################CV                                               
  CV<-numeric(dp.n)
  Wt<-matrix(numeric(dp.n*dp.n),ncol=dp.n)
  for (i in 1:dp.n)
  {
    if (DM.given)
         dist.vi<-dMat[,i]
    else
    {
       dist.vi<-gw.dist(dp.locat=dp.locat, focus=i, p=p, theta=theta, longlat=longlat)         
    }
    W.i<-gw.weight(dist.vi,bw,kernel,adaptive)
    #W.i<-gwr.Gauss(dist.vi^2, bw)
    #print(W.i)
    W.i[i]<-0
    Wt[,i]<-W.i
  }
  wt2 <- rep(1,dp.n)
  if (family=="poisson")
  {
     res1 <- gwr.poisson.wt(Y,X,bw,Wt)
     wt2<-res1[[1]]
     y.adj <- res1[[3]]
  } 
  else if (family=="binomial")
  {
     res1 <- gwr.binomial.wt(Y,X,bw,Wt)
     wt2<-res1[[1]]
     y.adj <- res1[[3]]   
  }
  for (i in 1:dp.n)
  {
    ##lm.i <- try(lm.wfit(y = y, x = x, w = w.i))  
    W.i<-Wt[,i]*wt2
    fun1<-function(X,Y,W.i) {betai<- solve(t(X*W.i)%*%X)%*%{t(X*W.i)%*%Y}}
    gw.resi<-try(fun1(X,y.adj,W.i))
  
    #gw.resi <- try(lm.wfit(y = Y, x = X, w = W.i))

    if(!inherits(gw.resi, "try-error"))
    {
      #b <- coefficients(gw.resi)
      yhat.noi<-X[i,]%*%gw.resi
      #CV[i] <- Y[i] - (t(b) %*% X[i,])
      if (family=="poisson")
        CV[i]<-Y[i]- exp(yhat.noi)
      else if (family=="binomial")
        CV[i]<-Y[i]-exp(yhat.noi)/(1+exp(yhat.noi))
      #CV[i]<-Y[i]-yhat.noi
      
    }
    else
    {
      CV[i]<-Inf
      break
    }
  }
  if (!any(is.infinite(CV)))
     CV.score<-t(CV) %*% CV   ### why squared errors are evaluated here? (TN)
  else
     {
        CV.score<-Inf
     }  
  if(adaptive)
    cat("Adaptive bandwidth:", bw, "CV score:", CV.score, "\n")
  else
    cat("Fixed bandwidth:", bw, "CV score:", CV.score, "\n")
  
  CV.score
}
 #Contribution of each observation to the score statistic used in cross-validation for ggwr
ggwr.cv.contrib<-function(bw, X, Y,family="poisson", kernel="bisquare",adaptive=F, dp.locat, p=2, theta=0, longlat=F,dMat)
{
   dp.n<-length(dp.locat[,1])
   #########Distance matrix is given or not

  if (is.null(dMat))
      DM.given<-F
  else
  {
    DM.given<-T
    dim.dMat<-dim(dMat)
    if (dim.dMat[1]!=dp.n||dim.dMat[2]!=dp.n)
    stop ("Dimensions of dMat are not correct")
  }
  ############################################CV
  CV<-numeric(dp.n)
  Wt<-matrix(numeric(dp.n*dp.n),ncol=dp.n)
  for (i in 1:dp.n)
  {
    if (DM.given)
         dist.vi<-dMat[,i]
    else
    {
       dist.vi<-gw.dist(dp.locat=dp.locat, focus=i, p=p, theta=theta, longlat=longlat)         
    }
    W.i<-gw.weight(dist.vi,bw,kernel,adaptive)
    #W.i<-gwr.Gauss(dist.vi^2, bw)
    #print(W.i)
    W.i[i]<-0
    Wt[,i]<-W.i
  }
  wt2 <- rep(1,dp.n)
  if (family=="poisson")
  {
     res1 <- gwr.poisson.wt(Y,X,bw,Wt, verbose=F)
     wt2<-res1[[1]]
     y.adj <- res1[[3]]
  } 
  else if (family=="binomial")
  {
     res1 <- gwr.binomial.wt(Y,X,bw,Wt, verbose=F)
     wt2<-res1[[1]]
     y.adj <- res1[[3]]   
  }
  for (i in 1:dp.n)
  {
    ##lm.i <- try(lm.wfit(y = y, x = x, w = w.i))  
    W.i<-Wt[,i]*wt2
    fun1<-function(X,Y,W.i) {betai<- solve(t(X*W.i)%*%X)%*%{t(X*W.i)%*%Y}}
    gw.resi<-try(fun1(X,y.adj,W.i))
  
    #gw.resi <- try(lm.wfit(y = Y, x = X, w = W.i))

    if(!inherits(gw.resi, "try-error"))
    {
      #b <- coefficients(gw.resi)
      yhat.noi<-X[i,]%*%gw.resi
      #CV[i] <- Y[i] - (t(b) %*% X[i,])
      if (family=="poisson")
        CV[i]<-Y[i]- exp(yhat.noi)
      else if (family=="binomial")
        CV[i]<-Y[i]-exp(yhat.noi)/(1+exp(yhat.noi))
      #CV[i]<-Y[i]-yhat.noi
      
    }
    else
    {
      CV[i]<-Inf
      break
    }
  } 
  CV
}
####Calculate the AICc with a given bandwidth
##Author: Binbin Lu
ggwr.aic<-function(bw, X, Y,family, kernel,adaptive, dp.locat, p=2, theta=0, longlat=F,dMat)
{
   dp.n<-length(dp.locat[,1])
   #########Distance matrix is given or not

  if (is.null(dMat))
      DM.given<-F
  else
  {
    DM.given<-T
    dim.dMat<-dim(dMat)
    if (dim.dMat[1]!=dp.n||dim.dMat[2]!=dp.n)
    stop ("Dimensions of dMat are not correct")
  }
  ############################################AIC
  ###In this function, the whole hatmatrix is not fully calculated and only the diagonal elements are computed
  S<-matrix(nrow=dp.n,ncol=dp.n)
  Wt<-matrix(numeric(dp.n*dp.n),ncol=dp.n)
  for (i in 1:dp.n)
  {
    if (DM.given)
         dist.vi<-dMat[,i]
    else
    {
       dist.vi<-gw.dist(dp.locat=dp.locat, focus=i, p=p, theta=theta, longlat=longlat)         
    }
    W.i<-gw.weight(dist.vi,bw,kernel,adaptive)
    Wt[,i]<-W.i
  }
  wt2 <- rep(1,dp.n)
  if (family=="poisson")
  {
     gw.possion.res<-gwr.poisson.wt(Y,X,bw,Wt)
     wt2<-gw.possion.res[[1]]
     llik<-gw.possion.res[[2]]
  }
  else if (family=="binomial")
  {
     gw.binomial.res<-gwr.binomial.wt(Y,X,bw,Wt)
     wt2<-gw.binomial.res[[1]]
     llik<-gw.binomial.res[[2]]
  }
  for (i in 1:dp.n)
  {
    #Ci=solve(t(X*W.i)%*%X)%*%{t(X*W.i)}
    W.i<-Wt[,i]*wt2
    fun2<-function(X,W.i) {Ci<-solve(t(X*W.i)%*%X)%*%{t(X*W.i)}}
    Ci<-try(fun2(X,W.i))
    #Ci<-solve(t(X*W.i)%*%X)%*%{t(X*W.i)}
   # gw.resi<-gw.reg(X,Y,W.i,hatmatrix=T,focus=i)
    #betas[i,]<-gw.resi[[1]] ######See function by IG
    #S[i,]<-gw.resi[[2]]
    if(!inherits(Ci, "try-error"))
      S[i,]<-X[i,]%*%Ci   
    else
    {
      S[i,]<-Inf
      break
    }  
  }
  
  if (!any(is.infinite(S)))
  {
    tr.S<-sum(diag(S))
    #AICc<--2*llik + 2*tr.S*dp.n/(dp.n-tr.S-2)
    AICc<--2*llik + 2*tr.S + 2*tr.S*(tr.S+1)/(dp.n-tr.S-1)   # This is generic form of AICc (TN)
  }
  else
    AICc<-Inf   
  if(adaptive)
    cat("Adaptive bandwidth (number of nearest neighbours):", bw, "AICc value:", AICc, "\n")
  else
    cat("Fixed bandwidth:", bw, "AICc value:", AICc, "\n")
  AICc
}


######### Two simplified fitting functions for bandwidth selection
############ Possipon GWGLM for bandwidth selection
gwr.poisson.wt<-function(y,x,bw,W.mat, verbose=T)
{
    ##Accuracy control
    tol<-1.0e-5
    maxiter<-20
    ############################################
	  var.n<-ncol(x)
    dp.n<-nrow(x)
    betas <-matrix(nrow=dp.n, ncol=var.n)
    betas1<- betas
    ##S: hatmatrix
    S<-matrix(nrow=dp.n,ncol=dp.n)
    ####################################
    ##model calibration
    it.count <- 0
    llik <- 0.0
    mu <- y + 0.1
    nu <- log(mu)
    if(verbose)
       cat(" Iteration    Log-Likelihood(With bandwidth: ",bw,")\n=========================\n")
    wt2 <- rep(1,dp.n)
    repeat {
     y.adj <- nu + (y - mu)/mu
     for (i in 1:dp.n)
     {
        W.i<-W.mat[,i]
        gw.resi<-gw.reg(x,y.adj,W.i*wt2,hatmatrix=F,i)
        betas1[i,]<-gw.resi[[1]]
     }
     nu <- gwr.fitted(x,betas1)
     mu <- exp(nu)
     old.llik <- llik
     llik <- sum(y*nu - mu - log(gamma(y+1)))
     if(verbose)
        cat(paste("   ",formatC(it.count,digits=4,width=4),"    ",formatC(llik,digits=4,width=7),"\n"))
     if (abs((old.llik - llik)/llik) < tol) break
     wt2 <- mu
     it.count <- it.count+1
     if (it.count == maxiter) break}
     res<-list(wt2,llik,y.adj)
     res    
}

############ Binomial GWGLM for bandwidth selection
gwr.binomial.wt<-function(y,x,bw,W.mat, verbose=T)
{
    tol=1.0e-5
    maxiter=20
    dp.n<-nrow(x)
    var.n <- ncol(x)
    betas <-matrix(nrow=dp.n, ncol=var.n)
    betas1<- betas
    ##S: hatmatrix
    S<-matrix(nrow=dp.n,ncol=dp.n) 
    ####################################
    ##model calibration
    n=rep(1,length(y))
    it.count <- 0
    llik <- 0.0
    mu <- 0.5
    nu <- 0
    if(verbose)
        cat(" Iteration    Log-Likelihood:(With bandwidth: ",bw,")\n=========================\n")
    wt2 <- rep(1,dp.n)
    repeat {
     y.adj <- nu + (y - n*mu)/(n*mu*(1 - mu))
     for (i in 1:dp.n)
     {
        W.i<-W.mat[,i]
        gw.resi<-gw.reg(x,y.adj,W.i*wt2,hatmatrix=F,i)
        betas1[i,]<-gw.resi[[1]]
     }
     nu <- gwr.fitted(x,betas1)
     mu <- exp(nu)/(1 + exp(nu))
     old.llik <- llik
     llik <- sum(lchoose(n,y) + (n-y)*log(1 - mu/dp.n) + y*log(mu/n))
     if(is.na(llik)) llik <-old.llik
     if(verbose)
        cat(paste("   ",formatC(it.count,digits=4,width=4),"    ",formatC(llik,digits=4,width=7),"\n"))
     if (abs((old.llik - llik)/llik) < tol) break
     wt2 <- n*mu*(1-mu)
     it.count <- it.count+1
     if (it.count == maxiter) break}
     res<-list(wt2,llik,y.adj)
     res    
}