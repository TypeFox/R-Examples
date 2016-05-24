plotIdentifiableZone <- function(nt,nb,add=FALSE,legend=TRUE,title=TRUE){
	ps<-seq(0,1,by=0.002)
	n<-length(ps)
	na<-nt-nb
	Nt<-nt*ps+2*(1-ps)
	fla<-na*ps/Nt    ## low limit of major allele, early
	fha<-(na*ps+1-ps)/Nt ## high limit major allele, early
	flb<-nb*ps/Nt    ## low limit of minor allele, early
	fhb<-(nb*ps+1-ps)/Nt ## high limit of minor allele, early
	flp<-rep(0,n)    ## low limit in aneuploid, late
	fhp<-ps/Nt   ## high limit in aneuploid, late
	flpp<-rep(0,n)   ## low limit in euploid
	fhpp<-(1-ps)/Nt  ## high limit in euploid
	col.list=rep(rgb(0,100,0,alpha=30,maxColorValue=255),4)
	if(!add){
        par(mar=c(4,4,2,0.1),xpd=TRUE)
        plot(0,0,cex=0,xlim=c(0,1),ylim=c(0,1),xlab='sAGP',ylab='SAF')
    }
	polygon(c(ps,ps[n:1]),c(fla,fha[n:1]),col=col.list[1],border=NA)
	polygon(c(ps,ps[n:1]),c(flb,fhb[n:1]),col=col.list[2],border=NA)
	polygon(c(ps,ps[n:1]),c(flp,fhp[n:1]),col=col.list[3],border=NA)
	polygon(c(ps,ps[n:1]),c(flpp,fhpp[n:1]),col=col.list[4],border=NA)
    col.borders=c('red2','orange3','dark green','deepskyblue')
	if(legend){
        legend(0.05,0.99,legend=c(expression(A[1]),expression(A[2]),'B','C'),lty=1,lwd=2,col=col.borders)
    }
    col.ccf=rgb(50,0,0,alpha=100,maxColorValue=255)
    if(nt==2&nb==0){
        polygon(c(0,0,1/3),c(0,0.5,1/3),col=col.ccf,border=NA)
    }
    delta=0.005
    if(nt==3&nb==1){
        xx=c(0,0,seq(0,0.5,by=0.002),seq(0.5,0,by=-0.002))
        yy=c(0,0.5,1/(2+seq(0,0.5,by=0.002)),2*seq(0.5,0,by=-0.002)/(2+seq(0.5,0,by=-0.002)))
        polygon(xx,yy,col=col.ccf,border=NA)
    }
    lines(ps+2*delta,fla+delta,col=col.borders[1],lwd=2)
    lines(ps+2*delta,fha+delta,col=col.borders[1],lwd=2)
    segments(2*delta,delta,2*delta,0.5+delta,col=col.borders[1],lwd=2)
    lines(ps+delta,flb+delta/2,col=col.borders[2],lwd=2)
    lines(ps+delta,fhb+delta/2,col=col.borders[2],lwd=2)
    segments(delta,delta/2,delta,0.5+delta/2,col=col.borders[2],lwd=2)
    lines(ps,flp-delta/2,col=col.borders[3],lwd=2)
    lines(ps,fhp-delta/2,col=col.borders[3],lwd=2)
    segments(1,0-delta/2,1,1/nt-delta/2,col=col.borders[3],lwd=2)
    lines(ps,flpp,col=col.borders[4],lwd=2)
    lines(ps,fhpp,col=col.borders[4],lwd=2)
    segments(0,0,0,0.5,col=col.borders[4],lwd=2)
    if(title)title(bquote(n[t]~'='~.(nt)~','~n[b]~'='~.(nb)),cex.main=2)
}