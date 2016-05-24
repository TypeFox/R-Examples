
svmmatch<-function(treat,X, burnin=100,gibbs=200,thin=2, dv=NULL){
	
	param=0#Vestigial variable
	if(is.null(dv)) dv<-as.matrix(rep(1,length(treat)))
			treat<-1*(treat==max(treat))
			y<-2*treat-1
			X<-X[,apply(X,2,sd)>0]
			X.orig<-X
			sds<-apply(X,2,FUN=function(x) sd(x[treat==1]))		
			X<-apply(X,2,FUN=function(x) (x-mean(x[treat==1]))/sd(x[treat==1]))
			svd1<-svd(X[treat==1,])
			X<-X%*%svd1$v
			sds2<-apply(X,2,FUN=function(x) sd(x[treat==1]))		
			X<-apply(X,2,FUN=function(x) (x-mean(x[treat==1]))/sd(x[treat==1]))

			b<-lm(y~X-1)$coef
			lambda<-as.vector(abs(1-y*(X%*%b))+.01)
		
		
			if(prod(dim(X))>1e06){
				n.keeps<-floor(1e06/ncol(X))
				keeps.cpp<-sample(1:nrow(X),n.keeps,replace=FALSE)
			} else{
				keeps.cpp<-(1:nrow(X))
				}
#bayesmatch_cpp <- cxxfunction(signature(X0="numeric",boldX0="numeric",lambda0="numeric",nu0="double",treat0="double", burnin0="int",total_gibbs0="int",thin0="int",param0="int",dv0="numeric"),   body=src1, plugin="RcppArmadillo")
		X0.cpp<-X[keeps.cpp,]
		boldX0.cpp<-X[keeps.cpp,]*y[keeps.cpp]*(treat[keeps.cpp]==0)	
		treat0.cpp<-	as.matrix(treat[keeps.cpp])
		dv0.cpp<-as.matrix(dv[keeps.cpp])
		lambda.cpp<-lambda[keeps.cpp]
		svm1<-bayesmatch_cpp(X0=X0.cpp,boldX0=boldX0.cpp, lambda0=lambda.cpp,nu0=(2/mean(b^2))^.5, treat0=treat0.cpp,burnin0=burnin,total_gibbs0=gibbs,thin0=thin,param,dv0=dv0.cpp)
		
		betas.rescale<-t(apply(svm1$b,1,FUN=function(x) x/sds2))
		betas.rescale<-betas.rescale%*%t(svd1$v)
		betas.rescale<-t(apply(svm1$b,1,FUN=function(x) x/sds))
		
		##Make balancing weights
		make.newwts<-function(b){
			fits<-X%*%b
			which.use<-(treat==1)|(1-y*fits>=0)
			mean1<-sum(which.use[treat==0&which.use&fits>0])
			mean0<-sum(which.use[treat==0&which.use&fits<0])

			wts2<-(which.use*1)
			wts2[fits>0&which.use==1]<-1/mean1
			wts2[fits<0&which.use==1]<-1/mean0
			wts2[treat==1&which.use==1]<-1

			wts2
		}

		wts2.all<-t(apply(svm1$b,1,make.newwts))
		effect.bal<-apply(wts2.all,1,FUN=function(x) lm(dv~treat,weights=x)$coef[2])
#system.time(effect.bal<-apply(wts2.all,1,FUN=function(x) lm(dv~treat,w=x)$coef[2]))

	

	
	return(
	list("effect"=effect.bal, "beta"=svm1$b_post, "margin"=svm1$marg, "bal.wts"=wts2.all,"X.scale"=X,"X.orig"=X.orig,"treat"=treat,"dv"=dv)
	)
}


#generate.simdata<-function(n,linear=FALSE){
#	k<-6
#	beta.Y1<-c(1,1,1,1,1,1)
#	beta.Y0<-c(1,-1,-1,-1,-1,1)
#	beta.T<-c(-2,2,1,2,1,2)
#
#	make.friedman<-function(x){
#		(10*sin(2*pi*x[,1]*x[,2])+20*(x[,3]-5)^2+10*x[,4]+5*x[,5])/200-4
#	}
#
#	X<-cbind(1,matrix(runif(n*(k-1)),nr=n))
#	expit<-function(x) (1+exp(-x))^-1
#
#	if(linear==TRUE){
#	Y1<-X%*%beta.Y1
#	Y0<- X%*%beta.Y0
#	}else{
#	Y1<-make.friedman(X[,-1])
#	Y0<--make.friedman(X[,-1])
#	}
#	#probs.true<-expit((X)%*%beta.T/5)
#	probs.true<-expit(make.friedman(X[,-1])*10-3)
#	treat.ind<-sample(1:length(probs.true), 100,pr = probs.true)
#	treat<-rep(0,length(probs.true))
#	treat[treat.ind]<-1
#
#	dv<-Y1*treat+Y0*(1-treat)
#	dv<-as.vector(dv)
#	lm(dv~treat)#Diff in Means approx .3
#	SATT<-mean(Y1[treat==1]-Y0[treat==1])
#	return(list("X"=X,"dv"=dv,"treat"=treat,"SATT"=SATT,"potential"=cbind(Y0,Y1),"propensity"=probs.true))
#}

balance<-function(treat,X,obj,plot.it=TRUE,sd.plot=.2,color=TRUE){
	
	par.old<-par()
	if(length(colnames(X))==0)  colnames(X)<-paste("X",1:ncol(X),sep="")
	X2<-apply(X,2,FUN=function(x) (x-mean(x[treat==1]))/sd(x[treat==1]))
	bal<-obj$bal.wts%*%X2
	bal.unwt<-colMeans(X2)
	if(plot.it){
		dens.all<-NULL
		for(i.dum in 1:ncol(X2)) {
			dens.all[[i.dum]]<-NULL
			dens.all[[i.dum]]<-density(bal[,i.dum])
		}
		col1<-rgb(0,0,1,alpha=.5)
		col2<-"red"
		col3<-"black"
		max.char<-max(sapply(colnames(X2),nchar))
		if(color==FALSE) col1<-col2<-col3<-"black"
		min.plot<-sapply(dens.all,FUN=function(obj) min(obj$x))
		max.plot<-sapply(dens.all,FUN=function(obj) max(obj$x))
		min.plot<-min(bal.unwt,min.plot)
		max.plot<-max(bal.unwt,max.plot)
		par(mar=c(3,max((max.char)/2.4,3),2,.3))
		plot(0,ylim=c(1,ncol(X2)+2),xlim=c(min.plot,max.plot),type="n",xlab="",ylab="",yaxt="n",xaxt="n")
		abline(h=seq(1,ncol(X2),1),col=gray(.75),lwd=1.25)
		segments(x0=0,x1=0,y0=0,y1=ncol(X2)+1.2,col=gray(.75))
		segments(x0=c(sd.plot,-sd.plot),x1=c(sd.plot,-sd.plot),y0=0,y1=ncol(X2)+1.2,lty=2,col=gray(.75))
		mtext(side=3,"Balance Plot")
		for(i in 1:ncol(X2)){
			lines(x=dens.all[[ncol(X2)-i+1]]$x,y=dens.all[[ncol(X2)-i+1]]$y/max(dens.all[[ncol(X2)-i+1]]$y)*.9+i,col=col1)
			text(bal.unwt[ncol(X2)-i+1],i+.3,"x",col=col3)
			text(colMeans(bal)[ncol(X2)-i+1],i,"|",font=2,cex=1,col=col2)
			text(quantile(bal[,ncol(X2)-i+1],c(0.025,.975)),i,":",font=2,col=col2)
			mtext(side=2,colnames(X2)[ncol(X2)-i+1],par("las"=1),at=(i+.5),cex=.8)
		}
		par("las"=0)
		axis(1)
		axis(1,at=c(sd.plot,-sd.plot))
		mtext(side=1,line=1.85,"Standardized Difference-in-Means",font=2)
#		abline(v=c(.1,-.1),lty=2)
#		abline(v=c(sd.plot,-sd.plot),lty=2)
		if(color==TRUE) col1<-rgb(0,0,1)
		legend("topleft",pch=c("x","|","|",":"),legend=c("Raw Data","Posterior Balance","Posterior Mean","95% CI"),horiz=TRUE,bty="n",col=c(col3,col1,col2,col2),x.intersp=.4)

	}
		par(par.old)
		output<-list("balance"=bal)
		invisible(output)
}


effect<-function(obj,color=TRUE,quant=c(0.025,0.975),legend.pos="topleft",label.main="Posterior Density of Effect Estimate",label.x="Outcome",label.y="Density"){

	par.old<-par()
	col1<-"red"
	col2<-rgb(0,0,1,alpha=.75)
	if(color==FALSE) col1<-col2<-FALSE
	plot(density(obj$effect),xlab="",ylab="",main="",col=col1,type="n")
	abline(h=0,col=gray(.7),lwd=1.5)
	abline(v=mean(obj$effect),lwd=2,col=col2)
	lines(density(obj$effect),col=col1,lwd=1.5)
	text(quantile(obj$effect,quant),0,"|",cex=1.5,col=col2)
	text(quantile(obj$effect,quant),0,"|",cex=1.5,col=col2)
	mtext(side=1,line=2,label.x,font=2)
	mtext(side=2,line=2,label.y,font=2)
	mtext(side=3,label.main,cex=1.25,line=.2,font=2)
	output<-list("quantiles"=quantile(obj$eff,seq(0,1,.05)),"mean"=mean(obj$eff))
	par(par.old)
	invisible(output)

}

autocorr<-function(obj){
	ac<-function(x) cor(x[-1],rev(rev(x)[-1]))
		apply(obj$beta,2,ac)
}

##Sensitivity analysis
sensitivity<-function(obj,seq.eval=seq(-1,1,.1),quant.eval=c(0.025,0.5,0.975),color=TRUE,legend.pos="topleft",
label.main="Sensitivity Analysis",label.x="Sensitivity Parameter",label.y="Outcome"){

par.old<-par()
fits.param<-obj$X.scale%*%t(obj$beta)
treat<-obj$treat
dv<-obj$dv
sens.wts<-function(fits,eps){
			y<-2*treat-1
			fits<-fits+y*eps
			which.use<-(treat==1)|(1-y*fits>=0)
			mean1<-sum(which.use[treat==0&which.use&fits>0])
			mean0<-sum(which.use[treat==0&which.use&fits<0])
			n<-length(treat)

			wts2<-(which.use*0)
			wts2[fits>0&which.use==1&treat==0]<- -1/mean(which.use&treat==0&fits>0)*.5
			wts2[fits<0&which.use==1&treat==0]<- -1/mean(which.use&treat==0&fits<=0)*.5
			#wts2[treat==0&which.use==1]<-wts2[treat==0&which.use==1]/mean(wts2[treat==0&which.use==1])*n
			wts2[treat==1]<-1/mean(treat)
			
			wts2[is.na(wts2)]<-0

			mean(wts2*dv)
			#lm(re78~treat,w=abs(wts2))$coef[2]
	}

sens.param<-sapply(seq.eval, FUN=function(eps.use) apply(fits.param,2,sens.wts,eps=eps.use))
sens.out<-apply(sens.param,2,FUN=function(x) quantile(x,quant.eval))

plot(0,xlim=range(seq.eval),ylim=range(sens.out),type="n",xlab="",ylab="",xaxt="n",yaxt="n")
axis(1)
axis(2)
abline(h=0,col=gray(.75),lwd=1.5)
abline(v=0,col=gray(.75),lwd=1.5)
col1<-rgb(1,0,0)
if(color==FALSE) col1<-"black"
lines(seq.eval,sens.out[1,],lty=2,col=col1,lwd=1.2)
lines(seq.eval,sens.out[2,],col=col1,lwd=1.5)
lines(seq.eval,sens.out[3,],lty=2,col=col1,lwd=1.2)
legend(legend.pos,legend=paste(rev(quant.eval)*100, "%"),lty=c(2,1,2),lwd=c(1.2,1.5,1.2),bty="n",col="red")
mtext(side=2,label.y,font=2,line=2)
mtext(side=1,label.x,font=2,line=2)
mtext(side=3,label.main,font=2,line=.2,cex=1.2)

output<-list("sens.mat"=sens.param)
	par(par.old)
	invisible(output)}#Closes out sensitivity analysis

##Control overlap
control.overlap<-function(obj,color=TRUE,label.main="Assessing Control Overlap",
label.x="Size of Control Set",label.y="Mass"){
	par.old<-par()
	col1<-"red"
	col2<-"blue"
	if(color==FALSE) col1<-col2<-"black"
	counts<-rowSums(obj$margin[,obj$treat==0])
	tab1<-table(counts)
	tab2<-tab1/sum(tab1)
	count.order<-sort(unique(counts))
	plot(0,type="n",xlim=range(counts),ylim=range(tab2),xlab="",ylab="",xaxt="n",yaxt="n")
	for(i in 1:length(count.order)){
		segments(x0=count.order[i],x1=count.order[i],y0=0,y1=tab2[i],col=col2)
		points(count.order[i],tab2[i],pch=19,col=col1)
	}
	axis(1)
	axis(2)
	mtext(side=2,label.y,font=2,line=2)
	mtext(side=1,label.x,font=2,line=2)
	mtext(side=3,label.main,font=2,line=.2,cex=1.2)
	output<-list("counts"=counts)
	par(par.old)
	invisible(output)
}

treatment.overlap<-function(obj,color=TRUE,thresh=.95){
	par.old<-par()
	treat<-obj$treat
	fits.param<-obj$X.scale%*%t(obj$beta)
	fits.treat<-fits.param[treat==1,]
	X.treat<-obj$X.orig[treat==1,]
	no.overlap<-1*(rowMeans(fits.treat>1)>thresh)
	overlap.out<-ifelse(no.overlap,"No match","Match")
	glm1<-glm(no.overlap~X.treat,family="binomial")
	print(table(overlap.out))
	print(summary(glm1))
	num.rows<-ceiling(ncol(X.treat)/4)
	par(mfrow=c(4,num.rows),mar=c(2,2,2,.5))
	col1<-"red"
	col2<-rgb(1,0,0,.1)
	if(color==FALSE) {
		col1<-"black"
		col2<-NULL
		}
	for(i in 1:ncol(X.treat)) 
	boxplot(X.treat[,i]~overlap.out,main=colnames(X.treat)[i],border=col1,col=col2)
output<-list("no.overlap"=no.overlap,"logit"=glm1)
	par(par.old)
	invisible(output)
}#Closes out treatment overlap