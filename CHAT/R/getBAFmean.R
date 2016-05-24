getBAFmean <- function(bb){
	cat('.')
	if(length(bb)<10)return(NA)
	ff<-quantile(bb,na.rm=T,c(0.01,0.99))
	bb<-bb[bb>=ff[1]&bb<=ff[2]]
	if(length(bb)<=10)return(NA)
	x<-density(bb,na.rm=T)$x
	y<-density(bb,na.rm=T)$y
	if(1){
		s.approx<-c(0.06,0.06)
		mu.approx<-c(0.05,0.95)
		startList<-list(a=1/sqrt(2*pi)/s.approx[1],b=s.approx[1],c=mu.approx[1],d=1/sqrt(2*pi)/s.approx[2],e=s.approx[2],f=mu.approx[2])
		while(1){
			fit2<-try(nls(y~a/b*exp(-(x-c)^2/2/b^2)+d/e*exp(-(x-f)^2/2/e^2),start=startList,trace=F),TRUE)
			if(class(fit2)=='nls'){break}
			mu.approx<-mu.approx+c(0.05,-0.05)
			if(abs(mu.approx[1]-mu.approx[2])<=0.001){break}
			startList<-list(a=1/sqrt(2*pi)/s.approx[1],b=s.approx[1],c=mu.approx[1],d=1/sqrt(2*pi)/s.approx[2],e=s.approx[2],f=mu.approx[2])
		}
		fit1<-try(nls(y~1/sqrt(2*pi)/b*exp(-(x-a)^2/2/b^2),start=list(b=0.06,a=0.5)))
		if(class(fit2)=='nls'& class(fit1)!='nls'){
			para.fit<-summary(fit2)$parameters
			c<-para.fit[3,1]
			f<-para.fit[6,1]
			return(abs(f-c)/2)
		}
		if(class(fit2)!='nls' & class(fit1)=='nls')return(0)
		if(class(fit2)!='nls'&class(fit1)!='nls')return(NA)
		if(class(fit1)=='nls'&class(fit2)=='nls'){
			para.fit1<-summary(fit1)
			para.fit2<-summary(fit2)
			para.fit<-summary(fit2)$parameters
			c<-para.fit[3,1]
			f<-para.fit[6,1]
			if(abs(c-0.5)<0.08|abs(f-0.5)<=0.08){
				if(abs(c-0.5)<0.08&abs(f-0.5)<0.08)return(abs(f-c)/2)
				else{
					if(abs(c-0.5)+abs(f-0.5)>0.2)return(0) else return(abs(f-c)/2)
				}
			}
			if(para.fit1$sigma<3*para.fit2$sigma)return(0)
			if(max(c,f)>1|min(c,f)<0|c+f>1.2)return(NA)
			return(abs(f-c)/2)
		}
	}
}