rwhatbufCat2<-function(rast, sites,bufsizes,att=1,asList=FALSE){
	coords<-coordinates(sites)
	res.list<-rep(list(NA),length(bufsizes))
	names(res.list)<-bufsizes
	for(i in 1:length(bufsizes)) {res.list[[i]]<-rep(list(NA),nrow(coords));names(res.list[[i]])<-1:nrow(coords)}
	compteur<-0
	t0<-Sys.time()
	cat("Number of buffers to compute: ",length(bufsizes)*nrow(coords),"\n")
		for(i in 1:length(bufsizes)){
			for (j in 1:nrow(coords)) {
			compteur<-compteur+1
			if ((j == 3) & (bufsizes[i] == bufsizes[1])) {
				esTime<-(Sys.time()-t0)*(length(bufsizes)*nrow(coords)/2)
				cat(paste("Time at start: ",t0,"\nApproximate run time: ", format(esTime),"\nApproximate end time: ",t0+esTime,"\n"),sep="")
			}
			if (isTRUE(all.equal(compteur/100,trunc(compteur/100)))) cat(compteur,"\n")
			flush.console()
				buff<-polycirc(bufsizes[i],coords[j,])
				buffP<-Polygon(buff)
				minirast<-readGDALbbox(rast,buffP,silent=TRUE)
				fullgrid(minirast)<-FALSE
				idx<-inout(coordinates(minirast),buff)
				res.list[[i]][[j]]<-table(minirast@data[,att][idx],useNA="ifany")
			}
		}
		
if (!asList) {
	cat("Building data.frame...\n")
	res2<-rep(list(NA),length(res.list))
	for(i in 1:length(res.list)) {
	res2[[i]]<-unlist(sapply(res.list[[i]],names))
	}
	res2<-sort(unique(unlist(res2)))
	mymat<-matrix(0,ncol=length(res2),nrow=length(bufsizes)*nrow(coords))
	colnames(mymat)<-res2
	compteur<-0
	for(i in 1:length(res.list)) {
		for (j in 1:nrow(coords)) {
			compteur<-compteur+1
			idx<-match(names(res.list[[i]][[j]]),colnames(mymat))
			mymat[compteur,idx]<-res.list[[i]][[j]]
		}
	}
	rep(bufsizes,each=nrow(coords))
  res.list<-data.frame(buffer=rep(bufsizes,each=nrow(coords)),ID=rep(1:nrow(coords),length(bufsizes)),mymat,check.names=FALSE)
}		
		
cat("Finished !\n")
return(res.list)
}
