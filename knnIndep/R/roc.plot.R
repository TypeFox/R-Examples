roc.plot <-
function(pows,legend=NULL,cols = colorRampPalette(c("blue","gray"))(dim(pows)[3]),mains = c("Linear","Quadratic","Cubic","Sine:period 1/2","Sine: period 1/8","X^(1/4)","Circle","Step function","Torus")){
	sens = seq(0,1,by=.1)
	sapply(1:dim(pows)[2],function(typ){
				plot(0,type="n",xlim=c(0,1),ylim=c(0,1),ylab=expression(alpha),xlab=expression(1-beta),main=mains[typ])
				spe=pows[,typ,]
				sapply(1:ncol(spe),function(i){points(1-spe[,i],sens,type="b",col=cols[i])})
				if(!is.null(legend)){
					legend("bottomright",legend=legend[,typ],col=cols,bg="white",lwd=1)
				}
			})	
}
