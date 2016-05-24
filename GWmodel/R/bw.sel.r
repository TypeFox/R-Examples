###Basic function for bandwidth selction
##Author: Binbin Lu
bw.gwr<-function(formula, data, approach="CV",kernel="bisquare",adaptive=FALSE, p=2, theta=0, longlat=F,dMat)
{
    ##Data points{
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
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0L)

  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  y <- model.extract(mf, "response")
  x <- model.matrix(mt, mf)
  dp.n<-nrow(data)
  ####################################Joke
  if(dp.n>1500)
  {
    cat("Take a cup of tea and have a break, it will take a few minutes.\n")
    cat("          -----A kind suggestion from GWmodel development group\n")
   }  
  #################### Recommond to specify a distance matrix
  if (missing(dMat))
  {
      DM.given<-F
      if(dp.n + dp.n <= 10000)
      {
        dMat <- gw.dist(dp.locat=dp.locat, rp.locat=dp.locat, p=p, theta=theta, longlat=longlat)
        DM.given<-T
      }
  }
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
       bw <- gold(gwr.cv,lower,upper,adapt.bw=adaptive,x,y,kernel,adaptive, dp.locat, p, theta, longlat,dMat)
    else if(approach=="aic"||approach=="AIC"||approach=="AICc")
       bw<-gold(gwr.aic,lower,upper,adapt.bw=adaptive,x,y,kernel,adaptive, dp.locat, p, theta, longlat,dMat)    
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
gwr.cv<-function(bw, X, Y, kernel="bisquare",adaptive=FALSE, dp.locat, p=2, theta=0, longlat=F,dMat, verbose=T)
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
    ##lm.i <- try(lm.wfit(y = y, x = x, w = w.i))
    fun1<-function(X,Y,W.i) {betai<- solve(t(X*W.i)%*%X)%*%{t(X*W.i)%*%Y}}
    gw.resi<-try(fun1(X,Y,W.i))
  
    #gw.resi <- try(lm.wfit(y = Y, x = X, w = W.i))

    if(!inherits(gw.resi, "try-error"))
    {
      #b <- coefficients(gw.resi)
      yhat.noi<-X[i,]%*%gw.resi
      #CV[i] <- Y[i] - (t(b) %*% X[i,])
      CV[i]<-Y[i]-yhat.noi
      
    }
    else
    {
      CV[i]<-Inf
      break
    }
  }
  if (!any(is.infinite(CV)))
     CV.score<-t(CV) %*% CV
  else
     {
        CV.score<-Inf
     }
  if(verbose)
  {
    if(adaptive)
      cat("Adaptive bandwidth:", bw, "CV score:", CV.score, "\n")
    else
      cat("Fixed bandwidth:", bw, "CV score:", CV.score, "\n")
  }
  CV.score
}

gwr.cv.contrib<-function(bw, X, Y, kernel="bisquare",adaptive=FALSE, dp.locat, p=2, theta=0, longlat=F,dMat)
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
    ##lm.i <- try(lm.wfit(y = y, x = x, w = w.i))
    fun1<-function(X,Y,W.i) {betai<- solve(t(X*W.i)%*%X)%*%{t(X*W.i)%*%Y}}
    gw.resi<-try(fun1(X,Y,W.i))
  
    #gw.resi <- try(lm.wfit(y = Y, x = X, w = W.i))

    if(!inherits(gw.resi, "try-error"))
    {
      #b <- coefficients(gw.resi)
      yhat.noi<-X[i,]%*%gw.resi
      #CV[i] <- Y[i] - (t(b) %*% X[i,])
      CV[i]<-Y[i]-yhat.noi
      
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
gwr.aic<-function(bw, X, Y, kernel="bisquare",adaptive=FALSE, dp.locat, p=2, theta=0, longlat=F,dMat, verbose=T)
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
  for (i in 1:dp.n)
  {
    if (DM.given)
         dist.vi<-dMat[,i]
    else
    {
       dist.vi<-gw.dist(dp.locat=dp.locat, focus=i, p=p, theta=theta, longlat=longlat)
    }
    W.i<-gw.weight(dist.vi,bw,kernel,adaptive)
    
    #Ci=solve(t(X*W.i)%*%X)%*%{t(X*W.i)}
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
    RSS.gw<-t(Y)%*%t(diag(dp.n)-S)%*%(diag(dp.n)-S)%*%Y
    sigma.hat2 <- RSS.gw/dp.n
    AICc<-dp.n*log(sigma.hat2) + dp.n*log(2*pi) + dp.n *((dp.n + tr.S) / (dp.n - 2 - tr.S))
  }
  else
    AICc<-Inf
  if(verbose)
  {     
    if(adaptive)
      cat("Adaptive bandwidth (number of nearest neighbours):", bw, "AICc value:", AICc, "\n")
    else
      cat("Fixed bandwidth:", bw, "AICc value:", AICc, "\n")
  }
  AICc
}


#######################################################################
#################### Golden Section search ############################
#######################################################################
# Golden selection function 
###Author: created by MC, edited by BL
# from www.me.gatech.edu/~aferri/me2016/Golden_Section.ppt
gold<-function(fun,xL,xU,adapt.bw=F,...){
   eps=1e-4 # working value - can be changed
   R <- (sqrt(5)-1)/2  # R <- 0.61803398....
   iter <- 1
   d <- R*(xU-xL)
   if (adapt.bw)
   {
     x1 <- floor(xL+d)
     x2 <- round(xU-d)
   }
   else
   {
     x1 <- xL+d
     x2 <- xU-d
   }  
   f1 <- eval(fun(x1,...))
   f2 <- eval(fun(x2,...))
   d1<-f2-f1
# Establish initial value of xopt:
   if (f1 < f2)
      xopt  <-  x1
   else xopt  <-  x2
# start main loop
   ea <- 100
   while ((abs(d) > eps) && (abs(d1) > eps)) {
      d <- R*d
      if   (f1 < f2) {
         xL  <-  x2
         x2  <-  x1
         if (adapt.bw)         
           x1 <- round(xL+d)
         else
           x1 <- xL+d
         #x1  <-  xL + d
         f2  <-  f1
         f1 <- eval(fun(x1,...))
         }
      else {
         xU  <-  x1
         x1  <-  x2
         if (adapt.bw)         
           x2 <- floor(xU - d)
         else
           x2  <-  xU - d 
         f1  <-  f2
         f2 <- eval(fun(x2,...))
         }
   iter  <-  iter + 1
# Establish value of xopt after iteration:
   if    (f1 < f2)
     xopt  <-  x1
   else xopt  <-  x2
   d1<-f2-f1
   ##print(paste(iter,f1,f2,xopt))
   }
    xopt
}