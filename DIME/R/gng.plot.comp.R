gng.plot.comp <-
function(data,obj, new.plot=TRUE,legpos=NULL, xlim=NULL, ylim=NULL, xlab=NULL, ylab=NULL,
  main=NULL,lwd=NULL,...){
obs <- unlist(data);
sNorm <- sort(obj$mu,index.return=TRUE);
sobs <- sort(obs);
len <- length(obs);
n <- matrix(0,len,obj$K);
I1 <- (obs < (-obj$th1))+0;
I2 <- (obs > (obj$th2))+0;
obsI1 <- obs*I1;
obsI11 <- obsI1[obsI1!=0];
obsI2 <- obs*I2;
obsI21 <- obsI2[obsI2!=0];
e1 <- obj$pi[1]*dexp(-sort(obsI11+obj$th1),1/obj$beta[1]);
e2 <- obj$pi[obj$K+2]*dexp(sort(obsI21-obj$th2),1/obj$beta[2]);
for (k in sNorm$ix){
  n[,k] <- obj$pi[k+1]*dnorm(sobs,obj$mu[k],obj$sigma[k]);
}
if (is.null(xlim)){
  xlim=c(sobs[1],sobs[length(sobs)]);
}
if(is.null(ylim)){
  ylim=c(0,max(e1,e2,n));
}
if(is.null(ylab)){
  ylab = "Density";
}
if(is.null(xlab)){
  xlab="Difference";
}
if(is.null(lwd)){
  lwd=2;
}
if(new.plot==TRUE){
  plot(sort(obsI11),e1,type="l",lty=obj$K+1,ylim=ylim,col=2,
    xlim=xlim,xlab=xlab,ylab=ylab,lwd=lwd,...);
  }else{
  lines(sort(obsI11),e1,type="l",lty=obj$K+1,ylim=ylim,col=2,
      xlim=xlim,xlab=xlab,ylab=ylab,lwd=lwd);  
}
  
tit1=NULL;
cols=NULL;
ltys = NULL;  
for (k in sNorm$ix){
  lines(sobs,obj$pi[k+1]*
    dnorm(sobs,obj$mu[k],obj$sigma[k]),type="l",col=k+3,lty=k,lwd=lwd);
  tmp = paste(round(obj$pi[k+1],2),
    " N(",prettyNum(obj$mu[k],digit=2),",",prettyNum(obj$sigma[k],digit=2),")",sep="");
  if((k+1) %in% obj$diffPiIdx){
    tmp = paste(tmp,"*",sep="");
  }
    cols = c(cols,k+3);
    ltys = c(ltys,k);
    tit1 <- c(tit1,tmp); 
}
lines(sort(obsI21),e2,type="l",col=3,lty=obj$K+2,lwd=lwd);

tit2 = paste(round(obj$pi[c(1,obj$K+2)],2)," E(",prettyNum(obj$beta,digit=2),")*",sep="");
cols <- c(cols,2,3);
ltys <- c(ltys,obj$K+1,obj$K+2);

if (is.null(legpos)){
  xleg = xlim[1];
  yleg = max(e1,e2,n);
}else{
  xleg = legpos[1];
  yleg = legpos[2];
}

if(new.plot==TRUE){
    legendTi = c(tit1,tit2[1],tit2[2]);
    legend(xleg,yleg,legendTi,col=cols,
      lty=ltys,lwd=rep(lwd,obj$K+2));
    if(is.null(main)) main="GNG";
    title(main=main);
  }else{
   legendTi = c("GNG",tit1,tit2[1],tit2[2]);
    legend(xleg,yleg,legendTi,col=c(1,cols),
      lty=c(1,ltys),lwd=c(3,rep(lwd,obj$K+2)));
  }
}

