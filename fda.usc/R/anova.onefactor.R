###################################################
anova.onefactor=function(object,group,nboot=100,plot=FALSE,verbose=FALSE,...){
#	require(MASS)
  fdata<-object
	if (class(fdata)!="fdata") stop("class(fdata) is not fdata")
	if (class(group)!="factor") group=as.factor(group)
	nda=ncol(fdata)
	nk=length(levels(group))
	if (verbose) cat("Number of groups",nk,"\n")
	index<-group==levels(group)[1]
	media=func.mean(fdata[index])
  total<-func.mean(fdata)
  n<-nrow(fdata)
  color<-numeric(n)
  color[index]<-2
	for (i in 2:nk){
	  index<-group==levels(group)[i]
    media<-c(media,func.mean(fdata[index]))
    color[index]<-i+1
    }
	cov<-array(dim=c(nda,nda,nk))
	for (i in 1:nk){cov[,,i]<-cov(fdata[[1]][group==levels(group)[i],])}
	wmest<-wmestadis(media)
	if (plot){  
	  dev.new(width=21,height=7)   
    palette("default")  
	  par(mfrow=c(1,3)) 
    mycols <- adjustcolor(palette("default"), alpha.f = 0.25)
    opal <- palette(mycols)
    plot(fdata,col=color,lty=2,lwd=.7,main="Functional Mean",...)
    palette("default")
    lines(c(total,media),col=1:(nk+1),lty=1,lwd=2)	  
#	  plot(c(total,media),lty=1,lwd=2,col=1:(nk+1),main="Functional Mean",...)
	  legend("bottom",c("Total",levels(group)),lwd=2,col=1:(nk+1),bty="n")
	  #	plot(density(wm),xlim=c(0,max(wm,wmest)),main=paste(" Groups:",nk,"\n","P-value=",pvalor," Stat:",round(wmest,3)))
	} 
	if (verbose) cat("Statistic:",wmest,"\n")
	media=c(media,total)	
	rangedat<-range(fdata)
	wm<-bootanovafunc(group,media,cov,nboot=nboot,plot=plot,rangedat)
  pvalor<-length(wm[wm>=wmest])/length(wm)
  if (plot){ 
  lines(c(total,media),col=1:(nk+1),lty=1,lwd=2)	   
  plot(density(wm),xlim=c(0,max(wm,wmest)),
       main=paste("Density plot for bootstrap resamples","\n","N levels:",nk,", Stat:",round(wmest,3),", P-value=",pvalor))  
  abline(v=wmest,col=2,lty=2,lwd=2)
#  points(wmest,0,pch=("*"),lwd=3)
  rug(wm)
		}
	return(list("pvalue"=pvalor,"stat"=wmest,"wm"=wm))
}

wmestadis<-function(xm){
	nk=nrow(xm)
	D=metric.lp(xm)
	estadis=sum(D[lower.tri(D)])
	return(estadis)
}
###################################################
 
bootanovafunc<-function(indi,xm,cov,nboot,plot=FALSE,rangedat=NULL){
nmues<-length(indi)
icont=table(indi)
nk=length(icont)
if (plot) {
  plot(xm,lwd=2,main=paste("Bootstrap resamples n=",nboot),lty=1,col=c(2:(nk+1),1),ylim=rangedat)
  legend("bottom",c("Total",levels(indi),"Resamples"),lwd=2,col=c(1,2:(nk+1),"grey"),bty="n")
  }
wm<-numeric(nboot)
xm3<-NULL
for(inb in 1:nboot){
	muestra=rproc2fdata(n=icont[1],t=argvals(xm),mu=xm[nk+1][[1]],sigma=cov[,,1])
	xm2=func.mean(muestra)
	for (k in 2:nk){
		muestra=rproc2fdata(n=icont[k],t=argvals(xm),mu=xm[nk+1][[1]],sigma=cov[,,k])
		aux=func.mean(muestra)
		xm2=c(xm2,aux)
			}
  if (inb>1) xm3<-c(xm3,xm2)
  else xm3<-xm2
if (plot) {lines(xm2,lty=2,col="grey",lwd=1)}
wm[inb]<-wmestadis(xm2)
}
return(wm)
}
###################################################

