image.sparseBC <-
function(x,labelx=TRUE,labely=TRUE,arrangex=FALSE,arrangey=FALSE,...){
res <- x$mus
dev.new(width=6,height=6)

	tempx <- NULL
	for(i in 1:max(x$Ds)){
		tempx <- c(tempx,which(x$Ds==i))
	}
	tempy <- NULL
	for(i in 1:max(x$Cs)){
		tempy <- c(tempy,which(x$Cs==i))
	}
	
	
	
if(arrangex==FALSE && arrangey==FALSE){
image(x=1:ncol(res),y=1:nrow(res),t(res),axes=TRUE,xaxt='n',yaxt='n',main="",xlab="",ylab="",col=heat.colors(12,alpha=0.8),...)
}

if(arrangex==FALSE && arrangey==TRUE){
	
image(x=1:ncol(res),y=1:nrow(res),t(res[tempy,]),axes=TRUE,xaxt='n',yaxt='n',main="",xlab="",ylab="",col=heat.colors(12,alpha=0.8),...)
}

if(arrangex==TRUE && arrangey==FALSE){
	
image(x=1:ncol(res),y=1:nrow(res),t(res[,tempx]),axes=TRUE,xaxt='n',yaxt='n',main="",xlab="",ylab="",col=heat.colors(12,alpha=0.8),...)
}

if(arrangex==TRUE && arrangey==TRUE){
	
image(x=1:ncol(res),y=1:nrow(res),t(res[tempy,tempx]),axes=TRUE,xaxt='n',yaxt='n',main="",xlab="",ylab="",col=heat.colors(12,alpha=0.8),...)
}


if(labelx==TRUE && arrangex==FALSE){
	axis(1,at=1:ncol(res),cex.axis=1.2)
}
if(labely==TRUE && arrangey==FALSE){
	axis(2,at=1:nrow(res),labels=c(nrow(res):1),cex.axis=1.2)
}

if(labelx==TRUE && arrangex==TRUE){

	axis(1,at=1:ncol(res),labels=c(tempx),cex.axis=1.2)
}
if(labely==TRUE && arrangey==TRUE){

	axis(2,at=1:nrow(res),labels=c(tail(tempy,nrow(res))),cex.axis=1.2)
}


title(main="Estimated Mean Matrix",xlab="Features",ylab="Observations",cex.lab=1.3,cex.main=1.5)

invisible(x)
invisible(res)

}
