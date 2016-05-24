SortModes <-
function(a.best,c.val=0.25){
  
ncp <- ncol(a.best$S);
Ng <- nrow(a.best$S);

#SORTING CRITERION

# A) Computation of relative data power. Store values in vector of size H=ncp. Could use different criterion here.
# need squared entries
 Ssq <- a.best$S * a.best$S ; 
 Asq <- a.best$A * a.best$A ;
 Xsq <- a.best$X * a.best$X ;
 rdp <- rep(0, times=ncp);
 for ( k in 1:ncp ) {
  rdp[k]<- sum(Ssq[,k])*sum(Asq[k,])/sum(Xsq) ;
 }
 rdp.s <- sort(rdp, na.last=NA,decreasing=TRUE, index.return=TRUE);

# B) sorting with mixture of contrast and data variance (Liebermeister)
 JG <- rep(0, times=ncp);
 JA <- rep(0, times=ncp);
 # generate values from std. normal distribution
 nu <- rnorm(10000,0,1);
 G0 <- mean(log(cosh(nu)));
 
 for ( k in 1:ncp ){
   # compute contrast for mode using logcosh
   G1 <- mean(log(cosh(a.best$S[,k])));
   JG[k] <- abs(G1-G0);
   JA[k] <- sum(Asq[k,]);
 }   
   sumJG <- sum(JG) ; sumJA <- sum(JA) ;
   J <- JG*(c.val/sumJG)+ JA*(1-c.val)/sumJA ;
   J.s <- sort(J, na.last=NA,decreasing=TRUE, index.return=TRUE);

   return(list(a.best=a.best,rdp=rdp.s,lbm=J.s));

} 