profPlot <-
function(specialization,id,col=c("black","red")){
	#specialization output from ordi.breadthvX
	cols<-rep(col[1],length(specialization$distances[[id]]))
	cols[which(specialization$group.vectors[id,]==TRUE)]<-col[2]
	
	specdist<-specialization$distances[[id]]#[which(specialization$group.vectors[id,]=="YES")]
	cols<-cols[order(specdist)]
	specdisto<-specdist[order(specdist)]
	
	plot(1:length(cols),specdisto,col=cols,pch=19,xlab="",ylab="distance",las=1,main=specialization$species[id])
}
