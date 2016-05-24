plot.epidemic <-
function(x, ...){
  opar<-par(mfrow=c(1,1))
  par(mfrow=c(1,1))
  x.data<-as.vector(as.matrix(x$i.data))
  semanas<-length(x.data)
  i.epi<-x$optimum.map[4]
  f.epi<-x$optimum.map[5]
	matplot(1:semanas,x.data,type="l",xlab="Week",ylab="Rate",col="#808080",lty=c(1,1),xaxt="n")
	if (!is.null(rownames(x$i.data))){
	  	axis(1,at=1:semanas,labels=rownames(x$i.data),cex.axis=1)
	}else{
		axis(1,at=1:semanas,labels=as.character(1:semanas),cex.axis=1)
	}
	if (is.na(i.epi)){
    puntos<-x.data
    points(1:semanas,puntos,pch=19,type="p",col="#00C000",cex=1.5)
	}else{
    # pre
    puntos<-x.data
    puntos[i.epi:semanas]<-NA
    points(1:semanas,puntos,pch=19,type="p",col="#00C000",cex=1.5)
    # epi
    puntos<-x.data
    if (i.epi>1) puntos[1:(i.epi-1)]<-NA
    if (f.epi<semanas) puntos[(f.epi+1):semanas]<-NA
    points(1:semanas,puntos,pch=19,type="p",col="#800080",cex=1.5)
    # post
    puntos<-x.data
    puntos[1:f.epi]<-NA
    points(1:semanas,puntos,pch=19,type="p",col="#FFB401",cex=1.5)
  }    
      
  legend(semanas*0.70,max.fix.na(x.data)*0.99,legend=c("Crude rate","Pre-epi period","Epidemic","Post-epi period"),
    lty=c(1,1,1,1),
    lwd=c(1,1,1,1),
    col=c("#808080","#C0C0C0","#C0C0C0","#C0C0C0"),
    pch=c(NA,21,21,21),
    pt.bg=c(NA,"#00C000","#800080","#FFB401"),
    cex=1)
  par(opar)
}
