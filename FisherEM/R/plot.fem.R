plot.fem <- function(x,frame=0,crit=c(),...){
old.par <- par(no.readonly = TRUE)
Y = eval.parent(x$call[[2]])
if ((frame==0 || frame==1) & length(x$plot$K)>1){
	K = x$plot$K; models = x$plot$models
        if (length(crit)==0) crit = x$crit
	if (crit=='bic') val = x$plot$bic; if (crit=='aic') val = x$plot$aic; if (crit=='icl') val = x$plot$icl 
	if (length(models)>1){ matplot(K,val,type='l',xlab='K',ylab=crit,col=1:length(models),lty=1:length(models),main='Selection of the number of groups',...)
		xleg = min(K) #+ 2/3*(max(K)-min(K))
		yleg = min(val) + 1/4*(max(val)-min(val))
		legend(xleg,yleg,models,col=1:length(models),lty=1:length(models),ncol=round(length(models)/3),cex=0.8,...)
	}
	else plot(K,val,type='b',xlab='K',ylab=crit,col=1:length(models),lty=1:length(models),main='Selection of the number of groups',...)
	if (x$call[[1]]=='sfem'){
		if (crit=='bic') val = x$plot$l1$bic; if (crit=='aic') val = x$plot$l1$aic; if (crit=='icl') val = x$plot$l1$icl 
		plot(x$plot$l1$l1,val,type='b',xlab='l1 value',ylab=crit,main='Selection of the sparsity penalty',...)
	}
}
if (frame==0 || frame==2) if (x$call[[1]]!='sfem') plot(x$loglik,type='b',xlab='Iterations',ylab='Log-likelihood',main='Log-likelihood',col=2,pch=20,cex=0.5,...)
if (frame==0 || frame==3){
	Y = as.matrix(Y)
	p = ncol(Y)
	n = nrow(Y)
	if (ncol(x$U)>1){
	  xx = as.matrix(Y)%*%x$U[,1:2]
	  cls = x$cls
	  min1= round(min(xx[,1]),1)-1
	  max1= round(max(xx[,1],1))+1
	  min2= round(min(xx[,2]),1)-1
	  max2= round(max(xx[,2]),1)+1
    topX = topY = 0
    xhist = yhist = list()
    for (k in 1:max(cls)){
      xhist[[k]] = hist(xx[cls==k,1],breaks=seq(min1,max1,0.1), plot=FALSE)
      yhist[[k]] = hist(xx[cls==k,2],breaks=seq(min2,max2,0.1), plot=FALSE)
      topX = max(topX,xhist[[k]]$counts); topY = max(topY,yhist[[k]]$counts)
    }
	  nf <- layout(matrix(c(2,0,1,3),2,2,byrow=TRUE), c(3,1), c(1,3), TRUE) 
	  xrange <- c(min1,max1); yrange <- c(min2,max2)
	  
	  par(mar=c(3,3,1,1)) 
	  plot(xx,col=cls,pch=cls,xlim=xrange, ylim=yrange,cex=1.5)
	  par(mar=c(0,3,1,1)) 
	  for (k in 1:max(cls)) barplot(xhist[[k]]$counts,axes=FALSE,col=k,ylim=c(0,topX),space=0,add=(k>1)) 
	  par(mar=c(3,0,1,1)) 
	  for (k in 1:max(cls)) barplot(yhist[[k]]$counts,axes=FALSE,horiz=TRUE,col=k,xlim=c(0,topY),space=0,add=(k>1))
	 }
   else {
	      cat('Since K=2, the data and the discriminative subspace have been projected on the 2 first PCs','\n')
		      if (ncol(Y)<nrow(Y))  Z = -eigen(cov(Y),symmetric=T)$vectors[,1:2]  
		      else {
			    z = eigen(Y%*%t(Y),symmetric=T)
			    Z = matrix(NA,p,2)
			    for (i in 1:2) Z[,i] = t(Y)%*%z$vectors[,i]/sqrt(n*z$values[i]) 
		      }
			  MU	= colMeans(Y)
			  proj= matrix(NA,2,p)
			  axU = matrix(x$U,p,1)%*%matrix(x$U,1,p)%*%matrix(10, p, 1)
			  proj= axU+matrix(MU,p,1)
			  Yproj = Y%*%Z
			  u 	= matrix(proj,1,p)%*%Z # projection sur les 2 pc
			  ybar 	= matrix(MU,1,p) %*% Z # proj des moyennes des cls sur les 2 pc
			  plot(Yproj,col=x$cls+1,xlab='comp.1',ylab='comp.2',pch=x$cls+1)
			  pente=(u[1,2]-ybar[1,2])/(u[1,1]-ybar[1,1])
			  oo=u[1,2]-pente*u[1,1]
			  xb=(2*ybar[1,1]-sqrt(50^2/(pente^2+1)))/2
			  xa=(2*ybar[1,1]+sqrt(50^2/(pente^2+1)))/2
        #cat(xa,xb,oo+pente*c(xa,xb),'\n')
			  lines(c(xa,xb),oo+pente*c(xa,xb),col=1,type='l',lwd=2)
        title(main='Estimated discriminative subspace (projected onto PCA axes)')
	}
}
par(old.par)
}

