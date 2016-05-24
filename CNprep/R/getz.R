getz <-
function(logr,emfit,zgroup,times){ 
	if(is.null(emfit$z))gz<-matrix(ncol=1,data=rep(1,length(logr)))
	else gz<-predict(emfit,newdata=logr)$z
	isfin<-matrix(ncol=ncol(gz),nrow=nrow(gz),data=is.finite(gz))
	gz[!isfin]<-0	#just being honest: we don't know how to assign these
	gz<-matrix(ncol=ncol(gz),
		data=apply(gz,2,cumsum)[seq(from=times,to=nrow(gz),by=times),])
	cisfin<-matrix(ncol=ncol(isfin),
		data=apply(isfin,2,cumsum)[seq(from=times,to=nrow(isfin),by=times),])
	cisfin<-cisfin-
		rbind(matrix(nrow=1,data=rep(0,ncol(cisfin))),cisfin[-nrow(gz),,drop=F])
	gz<-
		(gz-rbind(matrix(nrow=1,data=rep(0,ncol(gz))),gz[-nrow(gz),,drop=F]))/times
	gz[cisfin<times]<-NA
	return(gz%*%t(zgroup))
}
