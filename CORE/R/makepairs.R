makepairs <-
function(events){
	if("chrom"%in%dimnames(events)[[2]]){
		events<-events[order(events[,"chrom"]),,drop=F]
		chromstarts<-tapply(X=events[,"start"],INDEX=events[,"chrom"],FUN=min)
		chromends<-tapply(X=events[,"end"],INDEX=events[,"chrom"],FUN=max)
	}
	else{
		chromstarts<-min(events[,"start"])
		chromends<-max(events[,"end"])
	}
  ustarts<-sort(unique(events[,"start"])) #unique starts
  uends<-sort(unique(events[,"end"]))	#unique ends
	#find for each unique start and end its chromosome (sch and ech)	
  z<-cbind(c(chromstarts,ustarts),c(rep(1,length(chromstarts)),
      rep(0,length(ustarts))),c(rep(0,length(chromstarts)),1:length(ustarts)))
  z<-z[order(z[,1]),]
  sch<-cumsum(z[,2])[z[,3]!=0][order(z[z[,3]!=0,3])]
  z<-cbind(c(chromstarts,uends),c(rep(1,length(chromstarts)),
      rep(0,length(uends))),c(rep(0,length(chromstarts)),1:length(uends)))
  z<-z[order(z[,1]),]
  ech<-cumsum(z[,2])[z[,3]!=0][order(z[z[,3]!=0,3])]
	#for each unique end find the leftmost unique start in its chromosome
  sfirst<-match(1:length(chromstarts),sch)[ech]
	#for each unique end find the last preceding unique start
  z<-cbind(c(ustarts,uends),c(rep(1,length(ustarts)),rep(0,length(uends))))
  z<-z[order(z[,1]),]
  slast<-cumsum(z[,2])[z[,2]==0]
  z<-cbind(rep(1:length(uends),times=(slast-sfirst+1)),
        rep(0,sum(slast-sfirst+1)))
  z[match(1:length(uends),z[,1]),2]<-match(1:length(uends),z[,1])
  mypairs<-matrix(ncol=3,dimnames=list(NULL,c("start","end","id")),
		data=c(ustarts[(1:nrow(z))-cummax(z[,2])+sfirst[z[,1]]],
		uends[z[,1]],1:length(uends[z[,1]])))
  return(mypairs)
}
