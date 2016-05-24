timedelay.univariate.lm <-
function(bspline.data,target,regulator,maxdelay=ncol(bspline.data)*0.25,
single.adj.r.squared=0.8,min.coef=0.25,max.coef=4) {
	y<-bspline.data[target,]
	x<-bspline.data[regulator,]
	ts.point<-as.numeric(colnames(bspline.data))
	min.aic<- 10000
  best.adj.r.squared<- -10000  
	for (delay in 0:maxdelay) {
        x.to.fit<-x[1:(length(y)-delay)]
				y.to.fit<-y[(delay+1):length(y)]
				
				fit<-lm(y~.-1,data=data.frame(y=y.to.fit,x.to.fit))
				fit.aic<-AIC(fit)
				fit.coef<-summary(fit)$coef[,1]
								
				if ((fit.aic<min.aic) & sum(abs(fit.coef)>=min.coef,abs(fit.coef)<=max.coef)==2*length(fit.coef)) { 
					min.aic<-fit.aic
					best.adj.r.squared<- summary(fit)$adj.r.squared
				}   		
	} 
	if(best.adj.r.squared<single.adj.r.squared) return(NULL) 
  min.aic	
}

