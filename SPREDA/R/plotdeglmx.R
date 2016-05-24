plotdeglmx <-
function(x, type){
  if(class(x)!="deglmx"){print("Error: the object should be in the class of deglmx.")}
  dfs=x$dat$dfs
  dat=x$dat$dat
  dat.e=x$dyn.data
  dyn.mat=x$dyn.mat
  beta.index=x$beta.index
  coef=x$fit$coef[beta.index] #coef[1]=time.b 
  coef=c(0, coef)
  
  dyncovnames=x$dyncovnames
  n.dyn=length(dyncovnames)
  
  df.tmp2=cumsum(dfs)+2
  df.tmp1=c(2, df.tmp2[-n.dyn])+1

  if(missing(type)){
    for(i in 1:n.dyn){
      spline.mat=as.matrix(dyn.mat[, df.tmp1[i]:df.tmp2[i]])
      name=dyncovnames[i]
      eff=spline.mat%*%coef[df.tmp1[i]:df.tmp2[i]]
      obs=dat.e[, name]
      plot(obs[order(obs)],eff[order(obs)]-coef[2]/n.dyn, type="l",main="Effect plot",
           xlab=name,ylab=paste("Effect of ", name, sep=''), lwd=2)
      abline(h=0,lty=2)
    }
    ids=as.matrix(unique(dat[,1]))
    nn=length(ids)
    for(i in 1:nn){
      idx=(dat[,1]==ids[i])
      dat.i=dat[idx,]
      yhat.i=x$fit$fitted[idx]
      plot(dat.i[,2],dat.i[,3],type="p",main=ids[i], xlab="TIME", ylab="Observation")
      lines(dat.i[,2],yhat.i,col=2,lwd=2)
    } 
    
  }else if(type==1){
    for(i in 1:n.dyn){
      spline.mat=as.matrix(dyn.mat[, df.tmp1[i]:df.tmp2[i]])
      name=dyncovnames[i]
      eff=spline.mat%*%coef[df.tmp1[i]:df.tmp2[i]]
      obs=dat.e[, name]
      plot(obs[order(obs)],eff[order(obs)]-coef[2]/n.dyn, type="l",main="Effect plot",
           xlab=name,ylab=paste("Effect of ", name, sep=''), lwd=2)
      abline(h=0,lty=2)
    }
  }else if(type==2){
    ids=as.matrix(unique(dat[,1]))
    nn=length(ids)
    for(i in 1:nn){
      idx=(dat[,1]==ids[i])
      dat.i=dat[idx,]
      yhat.i=x$fit$fitted[idx]
      plot(dat.i[,2],dat.i[,3],type="p",main=ids[i], xlab="TIME", ylab="Observation")
      lines(dat.i[,2],yhat.i,col=2,lwd=2)
    } 
  }
  
}
