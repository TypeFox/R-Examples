gDist <-
function(out){
	
	# out : output of gPCA function
	pval<-ifelse(out$p.val<(1/out$nperm),paste("<",(1/out$nperm),sep=""),paste("=",out$p.val,sep=""))
if (is.null(out$filt)) {
	plot(density(out$delta.p),main=paste("Distribution of delta values from permutations\n(delta=",round(out$delta,3),"; p-value",pval,")",sep=""),
	xlim=c(min(out$delta.p,out$delta),max(out$delta.p,out$delta)))
	abline(v=out$delta,col='red',lty=2)
} else {
	
	plot(density(out$delta.p),main=paste("Distribution of delta values from permutations\n(delta=",round(out$delta,3),"; p-value",pval,"; filt=",out$filt,")",sep=""),
	xlim=c(min(out$delta.p,out$delta),max(out$delta.p,out$delta)))
	abline(v=out$delta,col='red',lty=2)
	}
}
