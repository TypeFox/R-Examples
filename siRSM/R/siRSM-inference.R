surface.test <- function(object)
{
   cat('Testing curvature of quadratic surface CONDITIONAL on the index ...\n\n')
   dxyz <- data.frame(u=object$u, v=object$v, y=object$y)
   lm0 <- lm(y~u+v,data=dxyz)
   lm1 <- lm(y~u+v+I(u*v),data=dxyz)
   if (object$type=='interaction.only') {
    return(anova(lm0, lm1))
   }
   else{
    lm2 <- lm(y~u+v+I(u^2)+I(u*v)+I(v^2),data=dxyz)
    anv1 <- anova(lm0, lm1, lm2)
	anv2 <- anova(lm0, lm2)
    return(list(anv1,anv2))
   }
}

#################################################################################
# ci.surface: bootstrap inference for various features of the response surface,
# CONDITIONAL on the estimated single index ...
#################################################################################
ci.surface <- function(obj, B=500, use.parallel=TRUE)
{   
  # calculate yhat
  S <- obj$u  
  E <- obj$v
  if (obj$type=='interaction.only'){
    beta=c(obj$coef[1],obj$coef[2],obj$coef[3],0,obj$coef[4],0)
  }
  else{
    beta=obj$coef
  }
  yhat <- cbind(rep(1,length(S)),S,E,S^2,S*E,E^2)%*%beta
  res=obj$residual
  
  if (B<=2500) { # parallel processing not worth the trouble
    use.parallel=FALSE
  }
  if (!use.parallel){
    cat('drawing',B,'parametric bootstrap samples ...\n\n')
  }
  if (use.parallel){  # new parallel code using foreach, etc
	MAX=detectCores()-1
    cl <- makeCluster(MAX,type='SOCK') 
	cat('Using parallel processing ...\n')
	print(cl)
	cat('\n')
    registerDoSNOW(cl=cl)
	clusterExport(cl, 'surface.stats.main')
    rslt <- foreach(b=1:B, .combine=rbind, .verbose=FALSE) %dopar% {
	    id <- sample(1:length(res),length(res),replace=T)
		res1 <- res[id]
		y1 <- yhat+res1;
		if (obj$type=='interaction.only'){
 		  betahat <- as.vector(coef(lm(y1~S+E+I(S*E))))
		  betahat <- c(betahat[1],betahat[2],betahat[3],0,betahat[4],0)
		}
		else{
 		  betahat <- as.vector(coef(lm(y1~S+E+I(S^2)+I(S*E)+I(E^2))))
		}
		paramf <- surface.stats.main(betahat)
		return(data.frame(b0=betahat[1],b1=betahat[2],b2=betahat[3],b3=betahat[4],b4=betahat[5],b5=betahat[6],
		    u0=paramf$u0,v0=paramf$v0,
			p10=paramf$p10, p11=paramf$p11,
			p20=paramf$p20, p21=paramf$p21,
			ax.congr=paramf$ax.congr, 
			ax2.congr=paramf$ax2.congr,
			ax.incongr=paramf$ax.incongr, 
			ax2.incongr=paramf$ax2.incongr))
		}
    stopCluster(cl)
   }
   else{
       rslt <- foreach(b=1:B, .combine=rbind, .verbose=FALSE) %do% {
	    id <- sample(1:length(res),length(res),replace=T)
		res1 <- res[id]
		y1 <- yhat+res1;
		if (obj$type=='interaction.only'){
 		  betahat <- as.vector(coef(lm(y1~S+E+I(S*E))))
		  betahat <- c(betahat[1],betahat[2],betahat[3],0,betahat[4],0)
		}
		else{
 		  betahat <- as.vector(coef(lm(y1~S+E+I(S^2)+I(S*E)+I(E^2))))
		}
		paramf <- surface.stats.main(betahat)
		return(data.frame(b0=betahat[1],b1=betahat[2],b2=betahat[3],b3=betahat[4],b4=betahat[5],b5=betahat[6],
		    u0=paramf$u0,v0=paramf$v0,
			p10=paramf$p10, p11=paramf$p11,
			p20=paramf$p20, p21=paramf$p21,
			ax.congr=paramf$ax.congr, 
			ax2.congr=paramf$ax2.congr,
			ax.incongr=paramf$ax.incongr, 
			ax2.incongr=paramf$ax2.incongr))
		}
   }
  
  # plots
  old=par()$mfrow
  par(mfrow=c(2,2))
  hist(rslt$ax.congr,xlab='',ylab='',main=expression(a[X]^c),col='grey',border=F)  
  hist(rslt$ax2.congr,xlab='',ylab='',main=expression(a[X^2]^c),col='grey',border=F)  
  hist(rslt$ax.incongr,xlab='',ylab='',main=expression(a[X]^i),col='grey',border=F)  
  hist(rslt$ax2.incongr,xlab='',ylab='',main=expression(a[X^2]^i),col='grey',border=F)  
#  plot(rslt$u0,rslt$v0, xlab='u', ylab='v', main='Stationary Point',pch='.',col='blue')
#  hist(rslt$p11,xlab=expression(p[11]),main='1st Prin. Ax.\nSlope',col='grey',border=F)
#  hist(rslt$p21,xlab=expression(p[21]),main='2nd Prin. Ax.\nSlope',col='grey',border=F)
#  hist(rslt$b0,xlab='',ylab='',main=expression(hat(beta)[0]),col='grey',border=F)  
#  hist(rslt$b1,xlab='',ylab='',main=expression(hat(beta)[1]),col='grey',border=F)  
#  hist(rslt$b2,xlab='',ylab='',main=expression(hat(beta)[2]),col='grey',border=F)  
#  hist(rslt$b3,xlab='',ylab='',main=expression(hat(beta)[3]),col='grey',border=F)  
#  hist(rslt$b4,xlab='',ylab='',main=expression(hat(beta)[4]),col='grey',border=F)  
#  hist(rslt$b5,xlab='',ylab='',main=expression(hat(beta)[5]),col='grey',border=F)  
  par(mfrow=old)
  
  # return CI  
  stationary.point=rbind(
  as.vector(c(quantile(rslt$u0,prob=0.025), quantile(rslt$u0,prob=0.975), mean(rslt$u0), sd(rslt$u0))),
  as.vector(c(quantile(rslt$v0,prob=0.025), quantile(rslt$v0,prob=0.975), mean(rslt$v0), sd(rslt$u0)))
  )
  stationary.point=as.data.frame(stationary.point)
  dimnames(stationary.point)[[1]]=c('u0','v0')
  dimnames(stationary.point)[[2]]=c('lower2.5','upper97.5','mean','std.err')

  prin.ax.1=rbind(
  as.vector(c(quantile(rslt$p10,prob=0.025), quantile(rslt$p10,prob=0.975), mean(rslt$p10), sd(rslt$p10))),
  as.vector(c(quantile(rslt$p11,prob=0.025), quantile(rslt$p11,prob=0.975), mean(rslt$p11), sd(rslt$p11)))
  )
  prin.ax.1=as.data.frame(prin.ax.1)
  dimnames(prin.ax.1)[[1]]=c('p10','p11')
  dimnames(prin.ax.1)[[2]]=c('lower2.5','upper97.5','mean','std.err')
  
  prin.ax.2=rbind(
  as.vector(c(quantile(rslt$p20,prob=0.025), quantile(rslt$p20,prob=0.975), mean(rslt$p20), sd(rslt$p20))),
  as.vector(c(quantile(rslt$p21,prob=0.025), quantile(rslt$p21,prob=0.975), mean(rslt$p21), sd(rslt$p21)))
  )
  prin.ax.2=as.data.frame(prin.ax.2)
  dimnames(prin.ax.2)[[1]]=c('p20','p21')
  dimnames(prin.ax.2)[[2]]=c('lower2.5','upper97.5','mean','std.err')

  beta=rbind(
  as.vector(c(quantile(rslt$b0,prob=0.025), quantile(rslt$b0,prob=0.975), mean(rslt$b0), sd(rslt$b0))),
  as.vector(c(quantile(rslt$b1,prob=0.025), quantile(rslt$b1,prob=0.975), mean(rslt$b1), sd(rslt$b1))),
  as.vector(c(quantile(rslt$b2,prob=0.025), quantile(rslt$b2,prob=0.975), mean(rslt$b2), sd(rslt$b2))),
  as.vector(c(quantile(rslt$b3,prob=0.025), quantile(rslt$b3,prob=0.975), mean(rslt$b3), sd(rslt$b3))),
  as.vector(c(quantile(rslt$b4,prob=0.025), quantile(rslt$b4,prob=0.975), mean(rslt$b4), sd(rslt$b4))),
  as.vector(c(quantile(rslt$b5,prob=0.025), quantile(rslt$b5,prob=0.975), mean(rslt$b5), sd(rslt$b5)))
  )
  beta=as.data.frame(beta)
  dimnames(beta)[[1]]=c('b0','b1','b2','b3','b4','b5')
  dimnames(beta)[[2]]=c('lower2.5','upper97.5','mean','std.err')
  
  line.congr=rbind(
  as.vector(c(quantile(rslt$ax.congr,prob=0.025), quantile(rslt$ax.congr,prob=0.975), mean(rslt$ax.congr), sd(rslt$ax.congr))),
  as.vector(c(quantile(rslt$ax2.congr,prob=0.025), quantile(rslt$ax2.congr,prob=0.975), mean(rslt$ax2.congr), sd(rslt$ax2.congr)))
  )
  line.congr=as.data.frame(line.congr)
  dimnames(line.congr)[[1]]=c('ax','ax2')
  dimnames(line.congr)[[2]]=c('lower2.5','upper97.5','mean','std.err')

  line.incongr=rbind(
  as.vector(c(quantile(rslt$ax.incongr,prob=0.025), quantile(rslt$ax.incongr,prob=0.975), mean(rslt$ax.incongr), sd(rslt$ax.incongr))),
  as.vector(c(quantile(rslt$ax2.incongr,prob=0.025), quantile(rslt$ax2.incongr,prob=0.975), mean(rslt$ax2.incongr), sd(rslt$ax2.incongr)))
  )
  line.incongr=as.data.frame(line.incongr)
  dimnames(line.incongr)[[1]]=c('ax','ax2')
  dimnames(line.incongr)[[2]]=c('lower2.5','upper97.5','mean','std.err')

  return(list(stationary.point=round(stationary.point,4),
              prin.ax.1=round(prin.ax.1,4),
			  prin.ax.2=round(prin.ax.2,4),
			  beta=round(beta,4),
			  line.congr=round(line.congr,4),
			  line.incongr=round(line.incongr,4)))
}

#################################################################################
# ci.index: bootstrap inference for the single index
#################################################################################
ci.index <- function(y, U, V, B=100, use.parallel=TRUE, ...)
{
  if (ncol(U) != ncol(V)) {
    cat('Covariate dimensions do not match!\n')
  }
  else{
   X=U; Z=V
   if (use.parallel){  # new parallel code using foreach, etc
	MAX=detectCores()-1
    cl <- makeCluster(MAX,type='SOCK') 
	cat('Using parallel processing ...\n')
	print(cl)
	cat('\n')
    registerDoSNOW(cl=cl)
	clusterExport(cl, c('siRSM', 'siRSM.default', 'single.run', 'multi.run'))
    w <- foreach(b=1:B, .combine=rbind, .verbose=FALSE) %dopar% {
	  ind=sample(1:nrow(X),size=nrow(X),replace=TRUE);
      m=siRSM(y[ind],X[ind,],Z[ind,],use.parallel=FALSE,...);
	  return(as.vector(m$w));
    }
	stopCluster(cl)
   }
   else{           # old serial code	
    w=matrix(0,B,ncol(X))
    for (b in 1:B){
	  cat('running bootstrap iteration',b,'\n')
	  ind=sample(1:nrow(X),size=nrow(X),replace=TRUE)
      m=siRSM(y[ind],X[ind,],Z[ind,],use.parallel=FALSE,...)
	  w[b,]=as.vector(m$w)
    }
   }	
    # plots
    old=par()$mfrow
	a=ceiling(sqrt(ncol(X)))
    par(mfrow=c(a,a))
    for (j in 1:ncol(X)){
      hist(w[,j],xlab='',ylab='',main=paste('w',j,sep=''),col='grey',border=F)  
    }
    par(mfrow=old)

    # return CI	
    index=data.frame(	
	lower2.5=apply(w,2,quantile,prob=2.5/100),
 	upper97.5=apply(w,2,quantile,prob=97.5/100),
	mean=apply(w,2,mean),
	std.err=apply(w,2,sd))
	dimnames(index)[[1]]=paste('w',1:ncol(X),sep='')
  }	
  return(round(index,4))
}
