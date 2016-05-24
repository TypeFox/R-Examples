randomRescaledEventMoves <-
function(events,boundaries,countOnly=F){
	chrlengths<-boundaries[,"end"]-boundaries[,"start"]+1
	chrfrac<-(events[,"end"]-events[,"start"]+1)/
		chrlengths[match(events[,"chrom"],boundaries[,"chrom"])]
	chrstart<-cumsum(chrlengths)-chrlengths+1
	newchrom<-sample(unique(events[,"chrom"]),size=nrow(events),replace=T)
	newstart<-chrstart[newchrom]+
		floor((1-chrfrac)*runif(nrow(events))*
		chrlengths[match(newchrom,boundaries[,"chrom"])])
	if(countOnly)return()
	newend<-newstart+
		pmax(0,round(chrfrac*chrlengths[match(newchrom,boundaries[,"chrom"])])-1)
	return(matrix(ncol=3,data=c(newstart,newend,newchrom),dimnames=list(NULL,
		c("start","end","chrom"))))
}
