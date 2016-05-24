#This is the main function available to the user
#This creates the object with some routine checks
#This also gives us the point estimates and can be used 
#with some summary and plot functions

ATE<- function(Y, Ti, X, theta=0, 
               ATT = FALSE, verbose = FALSE, 
               max.iter = 100, tol = 1e-10, initial.values = NULL,
               backtrack = TRUE, backtrack.alpha = 0.3, backtrack.beta = 0.5, ...){
  bt.a<- backtrack.alpha
  bt.b<- backtrack.beta
  
  if(is.vector(X)){
    X<- matrix(X, ncol = 1, nrow = length(X))
    warning("Data matrix 'X' is a vector, will be treated as n x 1 matrix")
  }
  if(nrow(X)!=length(Y)){
    stop("Dimensions of covariates and response do not match")
  }
  
  J<- length(unique(Ti))
  
  if(!all(0:(J-1) == sort(unique(Ti)))){
    stop("The treatment levels must be labelled 0,1,2,...")
  }
  #  if(!is.function(rho) | !is.function(rho1) | !is.function(rho2) ){
  #    stop("rho, rho1, rho2 must be user defined functions.
  #Alternatively, users can use the CR family of functions included in the package.")
  #  }
  #  if(!is.function(FUNu)){
  #    stop("'FUNu' must be a user defined function." )
  #  }
  if(!is.numeric(theta)){
    stop("theta must be a real number")
  }
  if(J>2 & ATT){
    stop("Treatment effect on the treated cannot be calculated for multiple treatment groups.")
  }
  
  K<- ncol(X)+1
  if(is.null(initial.values)){
    if(ATT){
      initial.values<- numeric(K)
    }else{
      initial.values<- matrix(0, ncol = K, nrow = J )
    } 
  }
  
  if(ATT & !is.vector(initial.values)){
    stop("For ATT we only need one vector of initial values for newton raphson")
  }
  if(!ATT){
    if( !is.matrix(initial.values) ){
      stop("Initial values must be a matrix")
    }
    if(any(dim(initial.values) !=  c(J,K)) ){
      stop("Matrix of initial values must have dimensions J x K.")
    }
  }
  
  rho<-function(x){cr.rho(x,theta=theta)}
  rho1<-function(x){d.cr.rho(x,theta=theta)}
  rho2<-function(x){dd.cr.rho(x,theta=theta)}
  
  FUNu<-  function(x) c(1,x)
  
  #Now we determine which category of the problem we are in
  gp<- "simple"
  if(ATT) gp<- "ATT"
  if(J>2) gp<- "MT"
  if(gp == "simple"){
    ini1<- initial.values[1,]
    ini2<- initial.values[2,]
    est<- get.est.simple(ini1, ini2, X, Y, Ti, rho, rho1, rho2,
                         FUNu, max.iter,
                         tol, backtrack, bt.a, bt.b,
                         verbose = verbose, ...)
    
  }else if(gp == "ATT"){
    ini2<- initial.values
    est<- get.est.ATT(ini2, X, Y, Ti, rho, rho1, rho2,
                      FUNu, max.iter,
                      tol, backtrack, bt.a , bt.b ,
                      verbose = verbose, ...)
  }else if(gp == "MT"){
    est<- get.est.MT(initial.values, X, Y, Ti, rho, rho1, rho2,
                     FUNu, max.iter,
                     tol, backtrack, bt.a, bt.b,
                     verbose, ...)
  }
  
  res<- est
  res$X<- X
  res$Y<- Y
  res$Ti<- Ti
  res$rho<- rho
  res$rho1<- rho1
  res$rho2<- rho2
  res$theta<- theta
  res$FUNu<- FUNu
  res$gp<- gp
  res$J<- J
  res$K<- K
  res$vcov<- estimate_variance(res,...)
  res$call<- match.call()
  if(gp=="simple"){
    est<- c(res$Y1,res$Y0,res$tau)
    names(est)<- c("E[Y(1)]", "E[Y(0)]", "ATE")
    res$est<- est
    res$Y1<- NULL
    res$Y0<- NULL
    res$tau<- NULL
  }else if(gp == "ATT"){
    est<- c(res$Y1,res$Y0,res$tau)
    names(est)<- c("E[Y(1)|T=1]", "E[Y(0)|T=1]", "ATT")
    res$est<- est
    res$Y1<- NULL
    res$Y0<- NULL
    res$tau<- NULL
  }else{
    est<- res$Yj.hat
    names(est)<- paste("E[Y(",0:(J-1),")]",sep = "")
    res$est<- est
    res$Yj.hat<- NULL
  }
  
  class(res)<- "ATE"
  return(res)
}

print.ATE<- function(x, ...){
  object<- x
  if(object$gp == "simple"){
    cat("Call:\n")
    print(object$call)
    cat("\nThe analysis was completed for a simple study design with binary treatment.\n")
    cat("\nPoint Estimates:\n")
    print(object$est)
  }else if(object$gp == "ATT"){
    cat("Call:\n")
    print(object$call)
    cat("\nThe analysis was completed for a simple study design with binary treatment.\n")
    cat("\nPoint Estimates:\n")
    print(object$est)
  }else{
    cat("Call:\n")
    print(object$call)
    cat("\nThe analysis was completed for a simple study design with binary treatment.\n")
    cat("\nPoint Estimates:\n")
    print(object$est)
  }
}


summary.ATE<- function(object, ...){
  if(object$gp== "simple" || object$gp== "ATT"){
    var.tau<- object$vcov[1,1]+object$vcov[2,2]-
      2*object$vcov[1,2]
    se<- c(sqrt(diag(object$vcov)), sqrt(var.tau))
  }else{
    se<- sqrt(diag(object$vcov))
  }
  Ci.l<- object$est+se*qnorm(0.025)
  Ci.u<- object$est+se*qnorm(0.975)
  
  z.stat<- object$est/se
  p.values<- 2*pnorm(-abs(z.stat))
  coef<- cbind(Estimate = object$est, 
               StdErr = se,
               "95%.Lower" = Ci.l,
               "95%.Upper" = Ci.u,
               Z.value = z.stat,
               p.value = p.values)
  
  if(object$gp == "simple"){
    res<- list(call = object$call, Estimate = coef, vcov = object$vcov,
               Conv = object$conv, Weights.p= object$weights.p,
               weights.q = object$weights.q)
    
  }else if(object$gp == "ATT"){
    res<- list(call = object$call, Estimate = coef, vcov = object$cov ,
               Conv = object$conv, weights.q = object$weights.q)
    
  }else{
    res<- list(call = object$call, Estimate = coef, vcov = object$cov ,
               Conv = object$conv, weights = object$weights.mat)
  }
  class(res)<- "summary.ATE"
  return(res)
  
}

print.summary.ATE <- function(x, ...){
  cat("Call:\n")
  print(x$call)
  cat("\n")
  printCoefmat(x$Estimate, P.values = TRUE, has.Pvalue=TRUE)
}

#####################################################################


my.ecdf<- function(t,x, weights = NULL){
  if(is.null(weights)){
    sum(1*(x <= t))/length(x)
  }else{
    sum(1*(x<=t)*weights)
  }
}

plot.ATE<- function(x, ...){
  object<- x
  Ti<- object$Ti
  if(object$gp == "simple"){
    w.p<- object$weights.p[Ti==1]
    w.q<- object$weights.q[Ti==0]
    
    names<- colnames(object$X)
    if(is.null(names)){
      p<- ncol(object$X)
      names<- paste("X",1:p,sep = "")
    }
    x1<- as.matrix(object$X[Ti==1,])
    x0<- as.matrix(object$X[Ti==0,])
    for(i in 1:ncol(x1)){
      if(i==2) par(ask = TRUE)
      
      if(length(unique(object$X[,i])) == 2){
        Treatment<- x1[,i]
        Placebo<- x0[,i]
        plot(c(0.5,1,2,2.5), c(2,mean(Treatment), mean(Placebo),2), pch = 16, cex = 1.5, 
             ylim = c(0,1), ylab = "Mean of group", xlab = "",col = c("blue","red"),
             main = "Unweighted", xaxt = "n")
        axis(side = 1, at = c(1,2), labels = c("Treatment", "Control") )
        abline(h = mean(object$X[,i]), lty = 2)
        new_treat<- sum(w.p*Treatment)
        new_control<- sum(w.q*Placebo)
        plot(c(0.5,1,2,2.5), c(2, new_treat, new_control,2), pch = 16, cex = 1.5, 
             ylim = c(0,1), ylab = "Mean of group", xlab = "",col = c("blue","red"),
             main = "Weighted", xaxt = "n")
        axis(side = 1, at = c(1,2), labels = c("Treatment", "Control") )
        abline(h = mean(object$X[,i]), lty = 2)
        
      }else{
        rng<- range(c(x1[,i],x0[,i]))
        my.seq<- seq(rng[1],rng[2],length = 100)
        temp1<- sapply(my.seq, my.ecdf,x = x1[,i])
        temp0<- sapply(my.seq, my.ecdf,x = x0[,i])
        par(mfrow = c(1,2))
        plot(my.seq,temp1,cex = 0.4,pch = 16,type = "l",lty = 1,col = "red",
             xlab = names[i],ylab = "empirical CDF",main = "Unweighted empirical CDF")
        lines(my.seq, temp0, cex = 0.4, pch = 16, lty = 2,col = "blue")
        legend("bottomright", c("Treatment", "Control"), lty = c(1,2), col = c("red", "blue"))
        
        temp1<- sapply(my.seq, my.ecdf,x = x1[,i],weights = w.p)
        temp0<- sapply(my.seq, my.ecdf,x = x0[,i], weights = w.q)
        plot(my.seq,temp1,cex = 0.4,pch = 16,type = "l",lty = 1,col = "red",
             xlab = names[i],ylab = "emperical CDF",main = "Weighted emperical CDF")
        lines(my.seq, temp0, cex = 0.4, pch = 16, lty = 2,col = "blue")
      }
    }
    par(ask = FALSE)
    
  }else if(object$gp == "ATT"){
    #w.p<- object$weights.p[Ti==1]
    w.q<- object$weights.q[Ti==0]
    
    names<- colnames(object$X)
    if(is.null(names)){
      p<- ncol(object$X)
      names<- paste("X",1:p,sep = "")
    }
    x1<- object$X[Ti==1,]
    x0<- object$X[Ti==0,]
    for(i in 1:ncol(x1)){
      if(i==2) par(ask = TRUE)
      
      if(length(unique(object$X[,i])) == 2){
        Treatment<- x1[,i]
        Placebo<- x0[,i]
        plot(c(0.5,1,2,2.5), c(2,mean(Treatment), mean(Placebo),2), pch = 16, cex = 1.5, 
             ylim = c(0,1), ylab = "Mean of group", xlab = "",col = c("blue","red"),
             main = "Unweighted", xaxt = "n")
        axis(side = 1, at = c(1,2), labels = c("Treatment", "Control") )
        abline(h = mean(x1[,i]), lty = 2)
        
        new_control<- sum(w.q*Placebo)
        plot(c(0.5,1,2,2.5), c(2, mean(Treatment), new_control,2), pch = 16, cex = 1.5, 
             ylim = c(0,1), ylab = "Mean of group", xlab = "",col = c("blue","red"),
             main = "Weighted", xaxt = "n")
        axis(side = 1, at = c(1,2), labels = c("Treatment", "Control") )
        abline(h = mean(x1[,i]), lty = 2)
        
      }else{
        rng<- range(c(x1[,i],x0[,i]))
        my.seq<- seq(rng[1],rng[2],length = 100)
        temp1<- sapply(my.seq, my.ecdf,x = x1[,i])
        temp0<- sapply(my.seq, my.ecdf,x = x0[,i])
        par(mfrow = c(1,2))
        plot(my.seq,temp1,cex = 0.4,pch = 16,type = "l",lty = 1,col = "red",
             xlab = names[i],ylab = "emperical CDF",main = "Unweighted emperical CDF")
        lines(my.seq, temp0, cex = 0.4, pch = 16, lty = 2,col = "blue")
        legend("bottomright", c("Treatment", "Control"), lty = c(1,2), col = c("red", "blue"))
        
        temp1<- sapply(my.seq, my.ecdf,x = x1[,i])
        temp0<- sapply(my.seq, my.ecdf,x = x0[,i], weights = w.q)
        plot(my.seq,temp1,cex = 0.4,pch = 16,type = "l",lty = 1,col = "red",
             xlab = names[i],ylab = "emperical CDF",main = "Weighted emperical CDF")
        lines(my.seq, temp0, cex = 0.4, pch = 16, lty = 2,col = "blue") 
      }
      
    }
    par(ask = FALSE)
    
  }else{
    
    wgt <- object$weights.mat
    names<- colnames(object$X)
    if(is.null(names)){
      p<- ncol(object$X)
      names<- paste("X",1:p,sep = "")
    }
    p<- ncol(object$X)
    for(i in 1:p){
      if(i==2) par(ask = TRUE)
      
      if(length(unique(object$X[,i])) == 2){
        J<- object$J
        
        plot(c(0:(J-1)), c(mean(object$X[Ti==0, i]), rep(2,J-1)), pch = 16, cex = 1.5, 
             ylim = c(0,1), ylab = "Mean of group", xlab = "Treatment group",col = 1,
             main = "Unweighted", xaxt = "n", xlim = c(-0.5,J-1+0.5))
        axis(side = 1, at = 0:(J-1), labels = paste(0:(J-1)) )
        abline(h = mean(object$X[,i]), lty = 2)
        for(j in 1:(object$J-1)){
          points(j, mean(object$X[Ti==j,i]), pch = 16, cex = 1.5, col = j+1)
        }
        
        plot(c(0:(J-1)), c(sum(object$X[Ti==0, i]*wgt[1,Ti==0]), rep(2,J-1)), pch = 16, cex = 1.5, 
             ylim = c(0,1), ylab = "Mean of group", xlab = "Treatment group",col = 1,
             main = "Weighted", xaxt = "n",xlim = c(-0.5,J-1+0.5))
        axis(side = 1, at = 0:(J-1), labels = paste(0:(J-1)) )
        abline(h = mean(object$X[,i]), lty = 2)
        for(j in 1:(object$J-1)){
          points(j, sum(object$X[Ti==j,i]*wgt[j+1,Ti==j]), pch = 16, cex = 1.5, col = j+1)
        }
        
      }else{
        rng<- range(object$X[,i])
        my.seq<- seq(rng[1],rng[2],length = 100)
        temp0<- sapply(my.seq, my.ecdf,x = object$X[object$Ti==0,i])
        par(mfrow = c(1,2))
        plot(my.seq,temp0,cex = 0.4,pch = 16,type = "l",lty = 1,col = 1,
             xlab = names[i],ylab = "emperical CDF",main = "Unweighted emperical CDF")
        for(j in 1:(object$J-1)){
          temp<- sapply(my.seq, my.ecdf,x = object$X[object$Ti==j,i])
          lines(my.seq, temp, cex = 0.4, pch = 16, lty = j+1 , col = j+1)
        }
        legend("bottomright", paste("gp",0:(object$J-1)), lty = 1:object$J, col = 1:object$J )
        
        temp0<- sapply(my.seq, my.ecdf,x = object$X[object$Ti==0,i], weights = wgt[1,Ti==0])
        plot(my.seq,temp0,cex = 0.4,pch = 16,type = "l",lty = 1,col = 1,
             xlab = names[i],ylab = "emperical CDF",main = "Weighted emperical CDF")
        for(j in 1:(object$J-1)){
          temp<- sapply(my.seq, my.ecdf,x = object$X[object$Ti==j,i], weights = wgt[j+1,Ti==j])
          lines(my.seq, temp, cex = 0.4, pch = 16, lty = j+1 , col = j+1)
        }
        
      }

    }
    par(ask = FALSE)
    
  }
}

estimate_variance<- function(object, ...){
  
  K<-  object$K
  if(object$gp == "simple"){
    cov.mat<- get.cov.simple(object$X, object$Y, object$Ti, object$FUNu, object$rho, 
                             object$rho1, object$rho2, object, ...)
    covEY<- cov.mat[-(1:(2*K) ),-(1:(2*K))]
    fin<- covEY
    #fin[3,3]<- covEY[1,1]+covEY[2,2] - 2*covEY[1,2]
    
  }else if(object$gp == "ATT"){
    cov.mat<-  get.cov.ATT(object$X, object$Y, object$Ti, object$FUNu, object$rho, 
                           object$rho1, object$rho2, object, ...)
    
    covEY<- cov.mat[-(1:K),-(1:K)]
    covEY<- covEY[-3,-3]
    fin<- covEY
    
  }else{
    cov.mat<-  get.cov.MT(object$X, object$Y, object$Ti, object$FUNu, object$rho, 
                          object$rho1, object$rho2, object, ...)
    covEY<- cov.mat[-(1:(object$J*K)),-(1:(object$J*K))]
    fin<- covEY
  }
  return(fin)
}
