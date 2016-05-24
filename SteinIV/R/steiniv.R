##################################################
### Trace operator.
tr <- function(X){
  out <- sum(diag(X));
  return(out)
}#tr

##################################################
### OLS
ols.est <- function(y,X,SE=FALSE){
   n  <- length(y);
   k  <- dim(X)[2];
   out <- list();
   XX <- t(X)%*%X; 
   iXX<- solve(XX);
   B  <- (iXX)%*%(t(X)%*%y);
   if(SE==TRUE){
     k       <- length(B);
     res     <- y - X%*%B;
     sig2.hat<- (t(res)%*%res)/(n-k);
     var.hat <- sig2.hat[1,1]*iXX;
     ### Out.
     out$est <- B;
     out$se  <- sqrt(diag(var.hat));
     out$var <- var.hat;
   }else{
     out$est <- B;
   }#fi
return(out)
}#ols.est

##################################################
### TSLS
tsls.est <- function(y,X,Z,SE=FALSE){
   n  <- length(y);
   k  <- dim(X)[2];
   out <- list();
   iZZ <- solve(t(Z)%*%Z);
   ZiZZ<- Z%*%iZZ;
   ZXh <- t(Z)%*%X;
   Xh  <- ZiZZ%*%ZXh;
   XXh <- t(Xh)%*%Xh; 
   iXXh<- solve(XXh);
   B   <- (iXXh)%*%(t(Xh)%*%y);
   if(SE==TRUE){
     k       <- length(B);
     res     <- y - X%*%B;
     sig2.hat<- (t(res)%*%res)/(n-k);     
     var.hat <- sig2.hat[1,1]*iXXh;
     ### Out.
     out$est <- B;
     out$se  <- sqrt(diag(var.hat));
     out$var <- var.hat;
   }else{
     out$est <- B;
   }#fi
return(out)
}#tsls.est

##################################################
### JIVE
jive.est <- function(y,X,Z,SE=FALSE,n.bt=100){
   n  <- length(y);
   k  <- dim(X)[2];
   out <- list();
   B   <- jive.internal(y,X,Z);
   ### SE
   if(SE==TRUE){
     k       <- length(B);
     var.hat <- matrix(0,k,k);
     for(b in 1:n.bt){
       bt <- sample(1:n,replace=TRUE);
       y.bt <- y[bt];
       X.bt <- as.matrix(X[bt,]);
       Z.bt <- as.matrix(Z[bt,]);
       B.bt <- jive.internal(y.bt,X.bt,Z.bt);
       var.hat <- var.hat + ((B.bt-B)%*%t(B.bt-B))/(n.bt-k); 
     }#n.bt     
     ### Out.
     out$est <- B;
     out$se  <- sqrt(diag(var.hat));
     out$var <- var.hat;
   }else{
     out$est <- B;
   }#fi
return(out)
}#jive.est

##################################################
### JIVE (Internal).
jive.internal <- function(y,X,Z){
   n <- length(y);
   k <- dim(X)[2];
   Xj <- matrix(0,n,k);
   iZZ<- solve(t(Z)%*%Z);
   Ga.hat <- (iZZ)%*%(t(Z)%*%X);
   h <- vector("numeric",n);
   for(i in 1:n){
      h[i]   <- t(Z[i,])%*%iZZ%*%Z[i,];    
      Xj[i,] <- (t(Z[i,])%*%Ga.hat - h[i]*X[i,])/(1-h[i]);    
   }#n
   XXj <- t(Xj)%*%X;
   iXXj<- solve(XXj);   
   #iXXj<- ginv(XXj);   ### Use Generalized inverse.
   B   <- (iXXj)%*%(t(Xj)%*%y);
return(B)
}#jive.internal

##################################################
### SPS
sps.est <- function(y,X,Z,SE=FALSE,ALPHA=TRUE,REF="TSLS",n.bt=100,n.btj=10){
   n <- length(y);
   k <- dim(X)[2];
   ### EST:
   out <- sps.internal(y,X,Z,REF=REF,ALPHA=ALPHA,n.btj=n.btj);
   #Verify: print(out);print(ALPHA)
   
   ### SE
   if(SE==TRUE){
     B       <- out$est;
     k       <- length(B);
     var.hat <- matrix(0,k,k);
     for(b.ind in 1:n.bt){
       bt <- sample(1:n,replace=TRUE);
       y.bt <- y[bt];
       X.bt <- as.matrix(X[bt,]);
       Z.bt <- as.matrix(Z[bt,]);
       B.bt <- sps.internal(y.bt,X.bt,Z.bt,REF=REF,ALPHA=FALSE,n.btj=n.btj)$est;
       var.hat <- var.hat + ((B.bt-B)%*%t(B.bt-B))/(n.bt-k); 
     }#n.bt
     ### Output
     out$se  <- sqrt(diag(var.hat));
     out$var <- var.hat
   }#fi   
return(out)
}#sps.est

##################################################
### CLS (Internal).
sps.internal <- function(y,X,Z,REF="TSLS",ALPHA=FALSE,n.btj=10){
 n <- length(y);
 k <- dim(X)[2];
 out <- list();
 
 ############################
 if(REF=="TSLS"){   
   ### OLS:
   XX <- t(X)%*%X; 
   iXX <- solve(XX);
   B.ols <- (iXX)%*%(t(X)%*%y);
   r.ols <- y - X%*%B.ols;
   sig2.ols<- (t(r.ols)%*%r.ols)/(n-k);
   var.ols <- sig2.ols[1,1]*iXX;

   ### TSLS:
   iZZ <- solve(t(Z)%*%Z);
   ZiZZ<- Z%*%iZZ;
   ZXh <- t(Z)%*%X;
   Xh  <- ZiZZ%*%ZXh;
   XXh <- t(Xh)%*%Xh; 
   iXXh<- solve(XXh);
   B.2ls <- (iXXh)%*%(t(Xh)%*%y);
   r.2ls <- y - X%*%B.2ls;
   sig2.2ls<- (t(r.2ls)%*%r.2ls)/(n-k);
   var.2ls <- sig2.2ls[1,1]*iXXh;

   ### CLS/SPS:
   sig2.cls<- (t(r.ols)%*%r.2ls)/(n-k);
   mse.ols.hat <- var.ols + (B.ols - B.2ls)%*%t(B.ols - B.2ls);
   cov.cls <- sig2.cls[1,1]*(iXX%*%(t(X)%*%Xh)%*%iXXh);
   alpha      <- tr(var.2ls-cov.cls)/tr(mse.ols.hat-2*cov.cls+var.2ls);
   B       <- alpha*B.ols + (1-alpha)*B.2ls; 
 }#TSLS 
 ############################
 if(REF=="JIVE"){
   ### OLS:
   XX <- t(X)%*%X; 
   iXX <- solve(XX);
   B.ols <- (iXX)%*%(t(X)%*%y);
   r.ols <- y - X%*%B.ols;
   sig2.ols<- (t(r.ols)%*%r.ols)/(n-k);
   var.ols <- sig2.ols[1,1]*iXX;
   
   ### JIVE:
   B.jive <- jive.internal(y,X,Z);

   ### BLCS:
   var.jive <- 0;
   cov.bcls <- 0;
   for(b in 1:n.btj){
      bt <- sample(1:n,replace=TRUE);
      y.bt <- y[bt];
      X.bt <- as.matrix(X[bt,]);
      Z.bt <- as.matrix(Z[bt,]);        
      B.ols.bt <- ols.est(y.bt,X.bt)$est;
      B.jive.bt<- jive.est(y.bt,X.bt,Z.bt)$est;
      var.jive <- var.jive + (B.jive.bt-B.jive)%*%t(B.jive.bt-B.jive)/(n.btj-k);
      cov.bcls <- cov.bcls + (B.jive.bt-B.jive)%*%t(B.ols.bt-B.ols)/(n.btj-k);
   }#n.btj
   #print(var.jive);print(cov.bcls)
   mse.ols.hat<- var.ols + (B.ols - B.jive)%*%t(B.ols - B.jive);
   alpha         <- tr(var.jive-cov.bcls)/tr(mse.ols.hat-2*cov.bcls+var.jive);
   B          <- alpha*B.ols + (1-alpha)*B.jive;
 }#JIVE
 
 ### ALPHA:
 if(ALPHA==TRUE){
    out$est <- B
    out$alpha  <- alpha
 }else{
    out$est <- B
 }#ALPHA 
return(out)
}#sps.internal
