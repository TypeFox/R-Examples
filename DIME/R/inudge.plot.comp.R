inudge.plot.comp <-
function(data,obj, new.plot=TRUE,legpos=NULL,xlim=NULL, ylim=NULL, xlab=NULL, ylab=NULL,
  main=NULL,lwd=NULL,...){
obs <- unlist(data);
if (is.null(xlim)){
  xlim=range(obs);
}
sNorm <- sort(obj$mu,index.return=TRUE);
sobs <- sort(obs);
len <- length(obs);
n <- matrix(0,len,obj$K);
d <- obj$pi[1]*dunif(sobs,obj$a,obj$b);
for (k in sNorm$ix){
  n[,k] <- obj$pi[k+1]*dnorm(sobs,obj$mu[k],obj$sigma[k]);
}
if (is.null(xlim)){
  xlim=c(sobs[1],sobs[length(sobs)]);
}
if(is.null(ylim)){
  ylim=c(0,max(n,d))
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
if(is.null(main)){
  if(obj$name=="iNUDGE"){
    main="iNUDGE";
  }else{
    main="NUDGE";}
}
if(new.plot){
  plot(sobs,d,type="l",col=2,lty=obj$K+1,xlim=xlim,ylim=ylim,xlab=xlab,
    ylab=ylab,...);
  }else{
  lines(sobs,d,type="l",col=2,lty=obj$K+1,xlim=xlim,ylim=ylim,xlab=xlab,
    ylab=ylab,lwd=lwd);
}

ti1=NULL;
cols=NULL;
ltys = NULL; 
for (k in sNorm$ix){
  lines(sobs,obj$pi[k+1]*
    dnorm(sobs,obj$mu[k],obj$sigma[k]),type="l",col=k+2,lty=k,lwd=lwd);
  tmp = paste(round(obj$pi[k+1],2), " N(",prettyNum(obj$mu[k],digit=2),", ",
    prettyNum(obj$sigma[k],digit=2),")",sep="");
  if((k+1) %in% obj$diffPiIdx){
    tmp = paste(tmp,"*",sep="");
  }
    cols = c(cols,k+2);
    ltys = c(ltys,k);
    ti1 <- c(ti1,tmp); 
}

ti2 = paste(round(obj$pi[1],2)," U(",prettyNum(obj$a,digit=2),",",
  prettyNum(obj$b,digit=2),")*",sep="");  
cols <- c(cols,2);
ltys <- c(ltys,obj$K+1);

if (is.null(legpos)){
  xleg = xlim[1];
  yleg = (9/10)*(max(n,d));
}else{
  xleg = legpos[1];
  yleg = legpos[2];
}

if (new.plot){
    legendTi = c(ti1,ti2);
    legend(xleg,yleg,legendTi,col=cols,
      lty=ltys,lwd=rep(lwd,obj$K+2));
    title(main=main)
  }else{
  if(obj$name=="iNUDGE"){
    legendTi = c("iNUDGE",ti1,ti2);
    legend(xleg,yleg,legendTi,col=c(1,cols),
      lty=c(1,ltys),lwd=c(3,rep(lwd,obj$K+2)));
  } else{
    legendTi = c("NUDGE",ti1,ti2);
    legend(xleg,yleg,legendTi,col=c(1,cols),
      lty=c(1,ltys),lwd=c(3,rep(lwd,obj$K+2)));
  }
  }
}

