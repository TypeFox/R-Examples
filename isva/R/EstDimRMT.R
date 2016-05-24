`EstDimRMT` <-
function(data.m,plot=TRUE){
 ### standardise matrix
 M <- apply(data.m,2,function(X){ (X - mean(X))/sqrt(var(X))});
 
 sigma2 <- var(as.vector(M));
 Q <- nrow(data.m)/ncol(data.m);
 ns <- ncol(data.m);
 lambdaMAX <- sigma2*(1+1/Q + 2*sqrt(1/Q));
 lambdaMIN <- sigma2*(1+1/Q - 2*sqrt(1/Q));
 delta <- lambdaMAX - lambdaMIN;#  print(delta);

 roundN <- 3;
 step <- round(delta/ns,roundN);
 while(step==0){
    roundN <- roundN+1;
    step <- round(delta/ns,roundN);
 }
  

 lambda.v <- seq(lambdaMIN,lambdaMAX,by=step);
 dens.v <- vector();
 ii <- 1;
 for(i in lambda.v){
   dens.v[ii] <- (Q/(2*pi*sigma2))*sqrt( (lambdaMAX-i)*(i-lambdaMIN) )/i;
   ii <- ii+1;
 }
 ## theoretical density
 thdens.o <- list(min=lambdaMIN,max=lambdaMAX,step=step,lambda=lambda.v,dens=dens.v);
 C <- 1/nrow(M) * t(M) %*% M;
 eigen.o <- eigen(C,symmetric=TRUE);
 ## empirical density
 estdens.o <- density(eigen.o$values,from=min(eigen.o$values),to=max(eigen.o$values),cut=0);
 intdim <- length(which(eigen.o$values > thdens.o$max));
 evalues.v <- eigen.o$values;
 ## plot
 if(plot){
  minx <- min(min(thdens.o$lambda),min(evalues.v));
  maxx <- max(max(thdens.o$lambda),max(evalues.v));
  miny <- min(min(thdens.o$dens),min(estdens.o$y));
  maxy <- max(max(thdens.o$dens),max(estdens.o$y));
  pdf("RMTplot.pdf",width=4,height=4);
  plot(thdens.o$lambda,thdens.o$dens,xlim=c(0.5,maxx),ylim=c(miny,maxy),type="b",col="green",xlab="Folded Eigenvalues",ylab="density",lwd=1.25);
  i <- min(which(estdens.o$x > min(evalues.v)));
  f <- max(which(estdens.o$x < max(evalues.v)));
  points(x=estdens.o$x[i:f],y=estdens.o$y[i:f],type="b",col="red",cex=0.5);
  for(i in 1:intdim){
   abline(v=evalues.v[i],col="red",lwd=2);
  }
  dev.off();
 }
 
 return(list(cor=C,dim=intdim,estdens=estdens.o,thdens=thdens.o,evals=eigen.o$values));
}

