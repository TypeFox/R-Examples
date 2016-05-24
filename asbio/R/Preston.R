Preston.dist<-function(counts, start=.2, cex.octave = 1, cex.legend = 1, cex.pt = 1,...){

	log2.counts<-logb(counts,2)
	hm<-hist(round(log2.counts,0),breaks=seq(0,round(max(log2.counts),0)),plot=FALSE)
	n<-hm$counts
	n.0<-max(hm$counts)
	octave.mode<-which(n==max(n))
	R<-(seq(0:(length(hm$counts)-1))-octave.mode)+.5
	no.spp<-nls(n~n.0*(exp(-1*(a*R)^2)),start=list(a=start))
	a<-summary(no.spp)$parameters[1]
	x<-seq(-1+-1*max(abs(c(min(R),max(R)))),1+max(abs(c(min(R),max(R)))))
	plot(x,c(rep(0,length(x)-1),n.0),type="n",xaxt="n",xlab=expression(paste(Log[2]," abundance (octaves)")),ylab="Number of species", ...)
	
	func<-function(x){n.0*exp(-1*(a*x)^2)}
	curve(func,from=min(x),to=max(x),add=TRUE,lwd=1)
	points(R,n,cex=cex.pt)
	staxlab(1,seq(min(R)-.5,max(R)+.5),hm$breaks,line.spacing=0, cex = cex.octave)
	abline(v=min(R)-0.5,lty=2,col=1)
	legend("topright",lty=2,col=1,legend=c("Veil line"),bty="n", cex = cex.legend)
	eq<-paste("n = ",bquote(.(n.0)),"*exp(",bquote(.(-1*round(a,3))),"*R)^2")
	Est.no.of.spp<-(n.0*sqrt(pi))/abs(a)
	res<-list(Equation=eq,Est.no.of.spp=Est.no.of.spp,Pct.sampling.completed=length(counts)/Est.no.of.spp*100)
	res
	}
