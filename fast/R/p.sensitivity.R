`p.sensitivity` <-
function(sen,
	xval,
	legend, 
        legend.cex=0.5,
	range=1:NROW(sen),
	col=1:NROW(sen),
	lty=rep(1,NROW(sen)),
	smooth=rep(FALSE,NROW(range)),
	x.range=1:length(xval), 
	m.max=max(sen[range,],na.rm=TRUE),
        limits=rep(FALSE,NROW(range)),
        xlab="time",
	ylab="Sensitivity",
	...){

stopifnot(mode(sen)=="numeric")
stopifnot(dim(sen)[2]>=max(x.range))

xval2 = xval[x.range]
sen2 = sen[,x.range]

plot(xval2,seq(0,m.max,length.out=NCOL(sen2)), col = 1,t="n",xlab=xlab, ylab=ylab, ...)
for(i in range){
	t.data<-as.double(sen2[i,])
	if(smooth[which(range==i)]){
		t.data <- lowess(t.data, f=1/2000)
		t.data <- t.data$y
	}
	lines(xval2,t.data , col = col[i],t="l", lty=lty[i])
        abline(h=limits[i], col = col[i], lty=2)
}

if(!is.na(eval(legend[1]))){
   legend("topright", inset=0.05, legend = legend[range],col=col[range], lty=lty[range], cex=legend.cex)
}

}

