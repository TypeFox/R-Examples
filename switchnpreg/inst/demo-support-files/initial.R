#############################################################
## Functions to calculate initial values for the functions ##
#############################################################

CVfunction.initial.Bayes <- function(S,x,y,sigma2){
              XDiff <- outer(x,x,FUN='-')
              K <- (1/(sqrt(2*pi)*S))*exp(-1/2*((XDiff)^2/(S^2)))
              A <- K+diag(rep(1,length(y)))*sigma2
              Ainv <- chol2inv(chol(A))
              SM <- K%*%Ainv
              n <- length(y)
              fhat <- SM%*%y
              GCV <- (1/n)*sum(((y-fhat)/(1-diag(SM)))^2)
              return(GCV)
              }
              
CVfunction.initial.PL <- function(lambda, y, B, R){
    
    #w <- rep(1/sigma2,length(y))
    #w <-  w/(sum(w)/length(y))
    S <-  B %*%solve(t(B) %*% B + 2 * exp(lambda) * R) %*% t(B) 
    fhat <- S %*% y
    
    gcv_i <- ((y-fhat) / (1-diag(S)))^2
    GCV <- (1/length(y)) * sum(gcv_i)
    
    GCV
      
}
              
fhatFunction.star <- function(x,x.star,y,S,sigma2){
                  XDiff <- outer(x,x,FUN='-')
                  K <- (1/(sqrt(2*pi)*S))*exp(-1/2*((XDiff)^2/(S^2)))
                  XDiff.star <- outer(x,x.star,FUN='-')
                  K.star <- (1/(sqrt(2*pi)*S))*exp(-1/2*((XDiff.star)^2/(S^2)))
                  A <- K+diag(rep(1,length(y)))*sigma2
                  Ainv <- chol2inv(chol(A))
                  fhat <- t(K.star)%*%Ainv%*%y
                  return(fhat)}  
                  
fhatFunction <- function(x,y,S,sigma2){
                XDiff <- outer(x,x,FUN='-')
                K <- (1/(sqrt(2*pi)*S))*exp(-1/2*((XDiff)^2/(S^2)))
                A <- K+diag(rep(1,length(y)))*sigma2
                Ainv <- chol2inv(chol(A))
                f <- K%*%Ainv%*%y
                Tr <- sum(diag(K%*%Ainv))
                out <- list(fhat=f,Trace=Tr)
                return(out)}

initial.indep.Bayes <- function(x,y,Interval1,Interval2,sigma2,nSubInt) {

  meanSpline  <- smooth.spline(x,y,spar=0.8)$y

  z.tmp <- vector(mode='numeric',length(y))
  z.tmp[which(y <= meanSpline)] <- 1
  z.tmp[which(y > meanSpline)] <- 2
  
  Y1 <- y[which(z.tmp==1)]
  x1 <- x[which(z.tmp==1)]
  Y2 <- y[which(z.tmp==2)]
  x2 <- x[which(z.tmp==2)]

  xmin1 <- vector(mode='numeric',length=nSubInt)
  objective1 <- vector(mode='numeric',length=nSubInt)
  for(i in 1:nSubInt){
      limInf <- Interval1[1]+(Interval1[2]/nSubInt)*(i-1)
      limSup <- i* (Interval1[2]/nSubInt)
      subInterval <- c(limInf,limSup)
      optimum <- optimize(CVfunction.initial.Bayes,subInterval,maximum=FALSE,x=x1,y=Y1,sigma2=sigma2[1])
      xmin1[i] <- as.vector(optimum$minimum)
      objective1[i] <- as.vector(optimum$objective)
      }
  
  xmin2 <- vector(mode='numeric',length=nSubInt)
  objective2 <- vector(mode='numeric',length=nSubInt)
  for(i in 1:nSubInt){
    limInf <- Interval2[1]+(Interval2[2]/nSubInt)*(i-1)
    limSup <- i* (Interval2[2]/nSubInt)
    subInterval2 <- c(limInf,limSup)
    optimum <- optimize(CVfunction.initial.Bayes,subInterval2,maximum=FALSE,x=x2,y=Y2,sigma2=sigma2[2])
    xmin2[i] <- as.vector(optimum$minimum)
    objective2[i] <- as.numeric(optimum$objective)
    }
  
  bestS1 <- xmin1[which(objective1==min(objective1))]
  bestS2 <- xmin2[which(objective2==min(objective2))]

  if (bestS1 < 10) {bestS1 <- 15}
  if (bestS2 < 10) {bestS2 <- 15}
  
  f1.x1.all <- fhatFunction(x=x1,y=Y1,S=bestS1,sigma2=sigma2[1])
  f1.x1 <- f1.x1.all$fhat
  trace1 <- f1.x1.all$Trace
  f1.xstar <- fhatFunction.star(x=x1,x.star=x2,y=Y1,S=bestS1,sigma2=sigma2[1])
  f1.x <- vector(mode='numeric',length(x))
  f1.x[which(z.tmp==1)]  <- f1.x1
  f1.x[which(z.tmp==2)]  <- f1.xstar
  
  f2.x2.all <- fhatFunction(x=x2,y=Y2,S=bestS2,sigma2=sigma2[2])
  f2.x2 <- f2.x2.all$fhat
  trace2 <- f2.x2.all$Trace
  f2.xstar <- fhatFunction.star(x=x2,x.star=x1,y=Y2,S=bestS2,sigma2=sigma2[2])
  f2.x <- vector(mode='numeric',length(x))
  f2.x[which(z.tmp==2)]  <- f2.x2
  f2.x[which(z.tmp==1)]  <- f2.xstar
  
  sigma2_1 <- sum((f1.x1-Y1)^2)/(length(Y1)-trace1)
  
  sigma2_2 <- sum((f2.x2-Y2)^2)/(length(Y2)-trace2)
  
  pooled_sigma2 <- ((length(Y1)-trace1)*sigma2_1 + (length(Y2)-trace2)*sigma2_2)/(length(Y1)+length(Y2)- trace1 - trace2)

  out <- list(f=cbind(as.vector(f1.x),as.vector(f2.x)),lambda=c(bestS1,bestS2),
              sigma2 =  c(pooled_sigma2,pooled_sigma2))
  return(out)
  
  }   

  



  
initial.indep.PL <- function(x, y) {

    meanSpline  <- smooth.spline(x, y, spar=0.8)$y

    Y1 <- y[which(y <= meanSpline)]
    x1 <- x[which(y <= meanSpline)]
    Y2 <- y[which(y > meanSpline)]
    x2 <- x[which(y > meanSpline)]

    s1 <- smooth.spline(x1, Y1)
    s2 <- smooth.spline(x2, Y2)
    
    #s1 <- smooth.spline(x1, Y1,spar=0.8)
    #s2 <- smooth.spline(x2, Y2,spar=0.8)
    
    lambda1 <- s1$lambda
    lambda2 <- s2$lambda
    
    trace1 <- sum(s1$lev) ## lev contains the diagonal elements of the smoother (hat) matrix
    trace2 <- sum(s2$lev)
    
    sigma2_1 <- sum((s1$y-Y1)^2)/(length(Y1)-trace1)
  
    sigma2_2 <- sum((s2$y-Y2)^2)/(length(Y2)-trace2)
  
    pooled_sigma2 <- ((length(Y1)-trace1)*sigma2_1 + (length(Y2)-trace2)*sigma2_2)/(length(Y1)+length(Y2)- trace1 - trace2)
    
    lambda1 <- s1$lambda
    lambda2 <- s2$lambda
    
    fit1 <- predict(s1, x)$y
    fit2 <- predict(s2, x)$y

    out <- list(f=cbind(as.vector(fit1),as.vector(fit2)),lambda=c(lambda1,lambda2),
              sigma2 =  c(pooled_sigma2,pooled_sigma2))
    return(out)
}


initialz.kmeans <- function (NW,Y.tmp,X.tmp) {

z.initial <- vector(mode='numeric',length=length(X.tmp))

nw <- NW #initial number of sub-intervals
d <- ceiling(length(X.tmp)/nw)

meanSpline  <- smooth.spline(X.tmp,Y.tmp,spar=0.8)$y
res <- meanSpline - Y.tmp

first.W <- seq(from=1,to=1*d,by=1)
x.tmp <- X.tmp[first.W]
y.tmp <- res[first.W] 

set.seed(1)
Kmeans <- kmeans(y.tmp,centers=2)
z.tmp <- Kmeans$cluster
means.tmp <-  Kmeans$centers

Data.tmp <- as.data.frame(list(y=y.tmp,z=z.tmp))

result.tmp <- summary(aov(y~z,Data.tmp))
pvalue.tmp <- result.tmp[[1]]$'Pr(>F)'[1]

pvalueNotsignif <- 0
  if (pvalue.tmp <= 0.05) {
    if (means.tmp[1] <= means.tmp[2]){
    one.indices <- which(z.tmp==1)
    two.indices <- which(z.tmp==2)
    z.tmp[one.indices] <- 2
    z.tmp[two.indices] <- 1
    z.initial[first.W] <- z.tmp} else{z.initial[first.W] <- z.tmp}
  } else {pvalueNotsignif <- pvalueNotsignif + 1}


for (i in 2:nw){
      if (i < nw) {
      ithwindow <- ((i-1)*d+1):(d*i)
      x.tmp <- X.tmp[ithwindow]
      y.tmp <- res[ithwindow] } else { ithwindow <- ((i-1)*d+1):length(X.tmp)
      x.tmp <- X.tmp[ithwindow]
      y.tmp <- res[ithwindow] 
      }

Kmeans <- kmeans(y.tmp,centers=2)
z.tmp <- Kmeans$cluster
means.tmp <-  Kmeans$centers

Data.tmp <- as.data.frame(list(y=y.tmp,z=z.tmp))

result.tmp <- summary(aov(y~z,Data.tmp))
pvalue.tmp <- result.tmp[[1]]$'Pr(>F)'[1]

pvalueNotsignif <- 0
  if (pvalue.tmp <= 0.05) {
    if (means.tmp[1] <= means.tmp[2]){
    one.indices <- which(z.tmp==1)
    two.indices <- which(z.tmp==2)
    z.tmp[one.indices] <- 2
    z.tmp[two.indices] <- 1
    z.initial[ithwindow] <- z.tmp} else{z.initial[ithwindow] <- z.tmp}
  } else {pvalueNotsignif <- pvalueNotsignif + 1}
}

return(z.initial)
}


initial.Markov.PL.Kmeans <- function(x,y,interval,basis,R,K) {
     
    z.tmp <- initialz.kmeans(NW=K,Y.tmp=y,X.tmp=x)
  
    Y1 <- y[which(z.tmp==1)]
    x1 <- x[which(z.tmp==1)]
    Y2 <- y[which(z.tmp==2)]
    x2 <- x[which(z.tmp==2)]

    B1<- getbasismatrix(evalarg=x1, basisobj=basis, nderiv=0)
    B2<- getbasismatrix(evalarg=x2, basisobj=basis, nderiv=0)
    B <- getbasismatrix(evalarg=x, basisobj=basis, nderiv=0)
    
    lambda_1 <- optimize(CVfunction.initial.PL,
                            interval,
                            maximum=FALSE,
                            y=Y1,
                            B=B1,R=R)$minimum
    
    lambda_2 <- optimize(CVfunction.initial.PL,
                            interval,
                            maximum=FALSE,
                            y=Y2,
                            B=B2,R=R)$minimum

     
    S1 <- B1 %*%solve(t(B1)  %*% B1 + 2 * exp(lambda_1) * R) %*% t(B1) 
    fhat1 <- B%*%solve(t(B1) %*% B1 + 2 * exp(lambda_1) * R) %*% t(B1)%*%Y1
    
    
    S2 <-  B2 %*%solve(t(B2) %*% B2 + 2 * exp(lambda_2) * R) %*% t(B2)
    fhat2 <- B%*%solve(t(B2) %*% B2 + 2 * exp(lambda_2) * R) %*% t(B2) %*%Y2

    trace1 <- sum(diag(S1)) 
    trace2 <- sum(diag(S2))
    
    sigma2_1 <- sum((S1%*%Y1-Y1)^2)/(length(Y1)-trace1)
  
    sigma2_2 <- sum((S2%*%Y2-Y2)^2)/(length(Y2)-trace2)
  
    pooled_sigma2 <- ((length(Y1)-trace1)*sigma2_1 + (length(Y2)-trace2)*sigma2_2)/(length(Y1)+length(Y2)- trace1 - trace2)
  
    out <- list(f=cbind(as.vector(fhat1),as.vector(fhat2)),lambda=c(lambda_1,lambda_2), 
                  sigma2 =  c(pooled_sigma2,pooled_sigma2))
    return(out)
}

initial.Markov.Bayes.Kmeans <- function(x,y,Interval1,Interval2,sigma2,nSubInt,K) {

  z.tmp <- initialz.kmeans(NW=K,Y.tmp=y,X.tmp=x)
  
  Y1 <- y[which(z.tmp==1)]
  x1 <- x[which(z.tmp==1)]
  Y2 <- y[which(z.tmp==2)]
  x2 <- x[which(z.tmp==2)]

  xmin1 <- vector(mode='numeric',length=nSubInt)
  objective1 <- vector(mode='numeric',length=nSubInt)
  for(i in 1:nSubInt){
      limInf <- Interval1[1]+(Interval1[2]/nSubInt)*(i-1)
      limSup <- i* (Interval1[2]/nSubInt)
      subInterval <- c(limInf,limSup)
      optimum <- optimize(CVfunction.initial.Bayes,subInterval,maximum=FALSE,x=x1,y=Y1,sigma2=sigma2[1])
      xmin1[i] <- as.vector(optimum$minimum)
      objective1[i] <- as.vector(optimum$objective)
      }
  
  xmin2 <- vector(mode='numeric',length=nSubInt)
  objective2 <- vector(mode='numeric',length=nSubInt)
  for(i in 1:nSubInt){
    limInf <- Interval2[1]+(Interval2[2]/nSubInt)*(i-1)
    limSup <- i* (Interval2[2]/nSubInt)
    subInterval2 <- c(limInf,limSup)
    optimum <- optimize(CVfunction.initial.Bayes,subInterval2,maximum=FALSE,x=x2,y=Y2,sigma2=sigma2[2])
    xmin2[i] <- as.vector(optimum$minimum)
    objective2[i] <- as.numeric(optimum$objective)
    }
  
  bestS1 <- xmin1[which(objective1==min(objective1))]
  bestS2 <- xmin2[which(objective2==min(objective2))]

  if (bestS1 < 10) {bestS1 <- 15}
  if (bestS2 < 10) {bestS2 <- 15}
  
  f1.x1.all <- fhatFunction(x=x1,y=Y1,S=bestS1,sigma2=sigma2[1])
  f1.x1 <- f1.x1.all$fhat
  trace1 <- f1.x1.all$Trace
  f1.xstar <- fhatFunction.star(x=x1,x.star=x2,y=Y1,S=bestS1,sigma2=sigma2[1])
  f1.x <- vector(mode='numeric',length(x))
  f1.x[which(z.tmp==1)]  <- f1.x1
  f1.x[which(z.tmp==2)]  <- f1.xstar
  
  f2.x2.all <- fhatFunction(x=x2,y=Y2,S=bestS2,sigma2=sigma2[2])
  f2.x2 <- f2.x2.all$fhat
  trace2 <- f2.x2.all$Trace
  f2.xstar <- fhatFunction.star(x=x2,x.star=x1,y=Y2,S=bestS2,sigma2=sigma2[2])
  f2.x <- vector(mode='numeric',length(x))
  f2.x[which(z.tmp==2)]  <- f2.x2
  f2.x[which(z.tmp==1)]  <- f2.xstar
  
  sigma2_1 <- sum((f1.x1-Y1)^2)/(length(Y1)-trace1)
  
  sigma2_2 <- sum((f2.x2-Y2)^2)/(length(Y2)-trace2)
  
  pooled_sigma2 <- ((length(Y1)-trace1)*sigma2_1 + (length(Y2)-trace2)*sigma2_2)/(length(Y1)+length(Y2)- trace1 - trace2)
  
  out <- list(f=cbind(as.vector(f1.x),as.vector(f2.x)),lambda=c(bestS1,bestS2),
              sigma2 =  c(pooled_sigma2,pooled_sigma2))
  return(out)
  
  }     







