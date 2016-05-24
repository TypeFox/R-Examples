# The Minkovski approach for Distance metric selection in a GWR calibration
# function for selecting an 'optimal' distance metric as to calibrate a given GWR model in a 'best fit' way
#Author: Binbin Lu
mink.approach <- function(formula, data, criterion="AIC", bw, bw.sel.approach = "AIC",adaptive=F, kernel="bisquare",
                          p.vals=seq(from=0.25, to=8, length.out=32), p.inf = T,
                          theta.vals = seq(from=0, to=0.5*pi, length.out=10), verbose=F, nlower = 10)
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
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0)

  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  y <- model.extract(mf, "response")
  x <- model.matrix(mt, mf)
  dp.n<-nrow(data)
  var.n <- ncol(x)
  longlat <- F
  #################Run through all the values of p and theta
  
  if(p.inf)
    p.vals <- c(p.vals, Inf) 
  np <- length(p.vals)
  ntheta <- length(theta.vals)
  diag.df <- matrix(nrow=np*ntheta, ncol=4)
  coefs.all <- array(0, dim=c(np,ntheta, dp.n, var.n))
  tag <- 0
  bw.given <- T
  if(missing(bw))
     bw.given <- F
  for (i in 1:np)
  {
    p=p.vals[i]
    cat("  ==========  Try Minkovski distance with p: ", p)
    cat("  ==========\n")
    for (j in 1:ntheta)
    {
      tag <- tag+1
      theta <- theta.vals[j]
      cat("  --The coordinate system is rotated by an angle", theta, "in radian\n")
      dMat <- gw.dist(dp.locat=dp.locat, p=p, theta=theta, longlat=longlat)
      diag.df[tag, 1] <- p
      diag.df[tag, 2] <- theta
      if(!bw.given)
        bw<-bw.gwr1(x, y, dp.locat,approach=bw.sel.approach,kernel=kernel,adaptive=adaptive,dMat,verbose=verbose, nlower = nlower)
      diag.df[tag, 3] <- bw
      if(criterion=="AIC" || criterion=="AICc")
      {
        res1 <- gwr.aic1(bw, x, y, kernel,adaptive, dp.locat, dMat)
        diag.df[tag, 4] <- res1[[1]]
        coefs.all[i,j,,] <- res1[[2]]
        if(adaptive)
          cat("    Adaptive bandwidth (number of nearest neighbours): ",bw, ", AICc value: ", diag.df[tag, 4], "\n") 
        else
          cat("    Fixed bandwidth: ",bw, ", AICc value: ", diag.df[tag, 4], "\n") 
      }
      else
      {
        res1 <- gwr.cv1(bw, x, y, kernel,adaptive, dp.locat, dMat)
        diag.df[tag, 4] <- res1[[1]]
        coefs.all[i,j,,] <- res1[[2]]
        if(adaptive)
          cat("    Adaptive bandwidth (number of nearest neighbours): ",bw, ", CV score: ", diag.df[tag, 4], "\n") 
        else
          cat("    Fixed bandwidth: ",bw, ", CV score: ", diag.df[tag, 4], "\n")  
      }
        #res1 <- gwr.aic1(bw, x, y, kernel,adaptive, dp.locat, dMat)
#        diag.df[tag, 4] <- res1[[1]]
#        coefs.all[i,j,,] <- res1[[2]]
#        if(adaptive)
#          cat("    Adaptive bandwidth (number of nearest neighbours): ",bw, ", AICc value: ", diag.df[tag, 4], "\n") 
#        else
#          cat("    Fixed bandwidth: ",bw, ", AICc value: ", diag.df[tag, 4], "\n")        
    }
  }
  #colnames(diag.df) <- c("p", "theta", "bandwidth", "AICc")
  if(criterion=="AIC" || criterion=="AICc")
    colnames(diag.df) <- c("p", "theta", "bandwidth", "AICc")
  else
    colnames(diag.df) <- c("p", "theta", "bandwidth", "CV")
  res <- list(diag.df, coefs.all)
}
##Bandwidth selection
bw.gwr1<-function(x, y, dp.locat,approach="AIC",kernel="gaussian",adaptive=FALSE,dMat, verbose=F, nlower = 10)
{
  dp.n <-  nrow(dp.locat)
  if(adaptive)
  {
    upper <- dp.n
    lower <- nlower
  }
  else
  {
      upper<-range(dMat)[2]
      lower<-upper/5000
  }
  ########################## Now the problem for the golden selection is too computationally heavy
    #Select the bandwidth by golden selection
    bw<-NA
    if(approach=="cv"||approach=="CV")
       bw <- gold(gwr.cv,lower,upper,adapt.bw=adaptive,x,y,kernel,adaptive, dp.locat,dMat=dMat,verbose=verbose)
    else if(approach=="aic"||approach=="AIC"||approach=="AICc")
       bw <- gold(gwr.aic,lower,upper,adapt.bw=adaptive,x,y,kernel,adaptive, dp.locat,dMat=dMat,verbose=verbose)
    bw
}

####Calculate the AICc with a given bandwidth
##Author: Binbin Lu
gwr.aic1<-function(bw, X, Y, kernel,adaptive,dp.locat, dMat)
{
  dp.n<-length(dp.locat[,1])
  S<-matrix(nrow=dp.n,ncol=dp.n)
  var.n <- ncol(X)
  betas <- matrix(nrow=dp.n,ncol=var.n) 
  for (i in 1:dp.n)
  {
    dist.vi<-dMat[,i]
    W.i<-gw.weight(dist.vi,bw,kernel,adaptive)
    #Ci=solve(t(X*W.i)%*%X)%*%{t(X*W.i)}
    fun2<-function(X,W.i) {Ci <- solve(t(X*W.i)%*%X)%*%{t(X*W.i)}}
    Ci<-try(fun2(X,W.i))
    if(!inherits(Ci, "try-error"))
    {
      S[i,] <- X[i,]%*%Ci
      betas[i,] <- Ci%*%Y 
      tr.S <- sum(diag(S))
      RSS.gw <- t(Y)%*%t(diag(dp.n)-S)%*%(diag(dp.n)-S)%*%Y
      sigma.hat2 <- RSS.gw/dp.n
      AICc <- dp.n*log(sigma.hat2) + dp.n*log(2*pi) + dp.n *((dp.n + tr.S) / (dp.n - 2 - tr.S))
    }  
    else
    {
      AICc <- Inf
      break
    }  
  }
  res <- list(AICc, betas)
}

gwr.cv1 <- function(bw, X, Y, kernel,adaptive, dp.locat,dMat)
{
   dp.n<-length(dp.locat[,1])
   var.n <- ncol(X)
   betas <- matrix(nrow=dp.n,ncol=var.n) 
  CV<-numeric(dp.n)
  for (i in 1:dp.n)
  {
    dist.vi<-dMat[,i]
    W.i<-gw.weight(dist.vi,bw,kernel,adaptive)
    ##lm.i <- try(lm.wfit(y = y, x = x, w = w.i))
    fun1<-function(X,Y,W.i) {betai<- solve(t(X*W.i)%*%X)%*%{t(X*W.i)%*%Y}}
    betas[i,] <- try(fun1(X,Y,W.i))
    #gw.resi <- try(lm.wfit(y = Y, x = X, w = W.i))
    W.i[i]<-0
    gw.resi<-try(fun1(X,Y,W.i))
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
  CV.score<-t(CV) %*% CV
  res<- list(CV.score, betas)
}

###Visualize the results of minkowski approach in a plotted matrix
mink.matrixview <- function(diag.df, znm=colnames(diag.df)[4], criterion="AIC")
{
  p.vals <- levels(as.factor(diag.df[,1]))
  #diag.df[!is.finite(diag.df[,1]),1]<- 9
  diag.df[,2] <- round(diag.df[,2], 2)
  theta.vals <- as.character(levels(as.factor(diag.df[,2])))
  np <- length(p.vals)
  ntheta <- length(theta.vals)
  x<- 1:np
  y<- 1:ntheta
  idx <- which(diag.df[,1] == 2)[1]
  if(criterion == "AIC" || criterion=="AICc")
  {
    aic.cv <-diag.df[,4]
    eu.val <-aic.cv[idx]
  }
  else
  {
    aic.cv <- log(diag.df[,4])
    eu.val <- log(aic.cv[idx])
    znm <- paste(paste("log(", znm, sep = ""), ")", sep="")
  }  
  
  plot(c(0,0), col="white", xlim = c(0, np+0.4), ylim = c(0, ntheta+0.4), axes=F, 
   xlab="Values of p", ylab = "Rotated anglein in radian", mar = c(0,0,0,0),
    main=paste(paste(znm, " value of the standard GWR calibration is :"), round(eu.val,2)))
   abline(v=0.5, lty=3, lwd=2,col="grey")
   abline(h=0.5, lty=3, lwd=2,col="grey")
   for (i in 1:np)
   {
      abline(v=i+0.5, lty=3, lwd=2,col="grey")
      text(x=i, y=0, p.vals[i],srt=90, font=2,)
   }
   for (j in 1:ntheta)
   {
      abline(h=j+0.5, lty=3, lwd=2,col="grey")
      text(x=-0.25, y=j, theta.vals[j], font=2)
   } 
   tag <-0
   dif <- round(aic.cv-eu.val, 2)
   for (i in 1:np)
     for(j in 1:ntheta)
     {
       tag <- tag+1
       col <- "black"
       if(dif[tag] >3)
         col <- "red"
       if (dif[tag] <= -3)
          col <- "green"
       text(x=x[i], y=y[j], as.character(dif[tag]),srt=90, cex=9/ntheta, col=col)
     }
   min.val <- min(dif)
   idx <- which.min(dif)
   #draw a rectangle at the optimu p and theta
   x.min <- which(p.vals==diag.df[idx,1])
   y.min <- which(theta.vals==diag.df[idx,2])
   rect(x.min-0.5, y.min-0.5, x.min+0.5, y.min+0.5, col = NA, border = NULL, lwd = 2)
   if(min.val < 0)
     cat("The maximum reduction of ", znm, "is: ", abs(min.val), "by the pair of values (p = ", diag.df[idx,1], "theta= ", diag.df[idx,2], ")\n")
   else
     cat("There is no reduction of ", znm, "from a non-Euclidean specified Minkovski distance\n")  
}