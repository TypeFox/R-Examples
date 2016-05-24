timedelay.lm <-
function(bspline.data, expr.data, target, regulator,   
maxdelay=ncol(bspline.data)*0.25,single.adj.r.squared=0.8,multiple.adj.r.squared=0.9,min.coef=0.25,max.coef=4,
output=FALSE,topdf=FALSE,xlab='time point',ylab='log ratio') {

	regulator.list<-regulator.aic.list<- NULL
	for (i in 1:length(regulator)) {
			regulator.i<-timedelay.univariate.lm(bspline.data,target,regulator[i],maxdelay=maxdelay,
      single.adj.r.squared=single.adj.r.squared,min.coef=min.coef,max.coef=max.coef)  	
			if (!is.null(regulator.i)) {
				regulator.list<-c(regulator.list,regulator[i])
				regulator.aic.list<- c(regulator.aic.list,regulator.i)
			}
	}
	names(regulator.aic.list)<- regulator.list
  regulator<-names(sort(regulator.aic.list,decreasing=FALSE))

	if (length(regulator)==0) return(NULL)

	y<-bspline.data[target,]
	min.aic<- 10000
	best.fit<- NULL

  regulator.delay.i<- regulator.delay<- NULL
  for (i in 1:length(regulator)) {
			regulator.delay.i<- c(regulator.delay,0); names(regulator.delay.i)<- c(names(regulator.delay),regulator[i])
			for (delay in 0:maxdelay) {
				regulator.delay.i[length(regulator.delay.i)]<-delay

				x.to.fit<-matrix(0,length(regulator.delay.i),length(y)-max(regulator.delay.i))
				rownames(x.to.fit)<-names(regulator.delay.i)
				for (j in 1:nrow(x.to.fit)) {
					x.to.fit[j,]<-bspline.data[rownames(x.to.fit)[j],(max(regulator.delay.i)-regulator.delay.i[j]+1):(length(y)-regulator.delay.i[j])]
				}
				y.to.fit<-y[(max(regulator.delay.i)+1):length(y)]

				fit<-lm(y~.-1,data=data.frame(y=y.to.fit,t(x.to.fit)))
				fit.aic<-AIC(fit)
				fit.coef<-summary(fit)$coef[,1]

				if ((fit.aic<min.aic) & sum(abs(fit.coef)>=min.coef,abs(fit.coef)<=max.coef)==2*length(fit.coef)
				& nrow(summary(fit)$coef) == length(names(regulator.delay.i))) {
					min.aic<-fit.aic
					regulator.delay<-regulator.delay.i
					best.fit<-fit
				}
			}
  }
       	
	if (summary(best.fit)$adj.r.squared<multiple.adj.r.squared) return(NULL)
	regulator.coef<- summary(best.fit)$coef[,1]   	
	
	ts.point<-as.numeric(colnames(bspline.data))
	if (output) {
		if (topdf) {
			pdf(paste(target,'.pdf',sep=''),width=8,height=6)
		} else {
			x11()
		}
		x.to.fit<-matrix(NA,length(regulator.delay),length(y)-max(regulator.delay))
		rownames(x.to.fit)<-regulator<- names(regulator.delay)		
		for (j in 1:nrow(x.to.fit)) {
			x.to.fit[j,]<-bspline.data[rownames(x.to.fit)[j],(max(regulator.delay)-regulator.delay[j]+1):(length(y)-regulator.delay[j])]*regulator.coef[j]
    }
		x<- bspline.data[rownames(x.to.fit),]
		pal<-rainbow(round(length(regulator.delay)*1.33))
		plot(ts.point[1:length(y)],y,type='l',col='black',ylim=c(min(y,x.to.fit,x),max(y,x.to.fit,x)),lwd=4,xlab=xlab,ylab=ylab)
		title(paste('adj.r.squared =',format(summary(best.fit)$adj.r.squared,digits=3)))
		
		lines(ts.point[(max(regulator.delay)+1):length(y)],predict(best.fit),lty='dashed',col='black',lwd=4)
		
		original.ts.point<- as.numeric(colnames(expr.data))
    points(original.ts.point, expr.data[target,],pch=4,col='black')

    for (j in 1:nrow(x.to.fit)) {
			lines(ts.point[1:length(y)],bspline.data[rownames(x.to.fit)[j],],col=pal[j],lwd=2)
			lines(ts.point[(max(regulator.delay)+1):length(y)],x.to.fit[j,],col=pal[j],lty='dashed',lwd=2)
			points(original.ts.point, expr.data[rownames(x.to.fit)[j],],pch=4,col=pal[j])
		}  		
		    
		legend.text<-c(target,paste('predicted',target))
		legend.col<-rep('black',2)
		legend.lty<-rep(c('solid','dashed'),length(regulator.delay)+1)
		legend.lwd<-c(2,2,rep(2,length(regulator.delay)*2))
		for (i in 1:length(regulator.delay)) {
			legend.text<-c(legend.text,paste(names(regulator.delay)[i],': ',format(regulator.coef[i],digits=3),sep=''),
      paste('delay: ',format(regulator.delay[i]*(ts.point[2]-ts.point[1]),digits=3)))
			legend.col<-c(legend.col,rep(pal[i],2))
		}
		legend('topleft',legend=legend.text,col=legend.col,lty=legend.lty,lwd=legend.lwd,cex=0.7,bty='n')
		
		if (topdf) dev.off()
	}
	
	return(list(delay=regulator.delay*(ts.point[2]-ts.point[1]),fit=best.fit))
}

