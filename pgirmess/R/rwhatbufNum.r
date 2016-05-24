rwhatbufNum<-function(rast, sites,bufsizes,att=1){
	coords<-coordinates(sites)
	mylist<-rep(list(NA),length(bufsizes))
	names(mylist)<-bufsizes
	compteur<-0
	bufcomp<-0
	cat("Number of buffers to compute: ",length(bufsizes)*nrow(sites@data),"\n")
	for(i in 1:length(mylist)){
		mylist[[i]]<-rep(list(NA),nrow(coords))
		names(mylist[[i]])<-1:nrow(coords)
	}
		t0<-Sys.time()
		for(i in bufsizes){
		bufcomp<-bufcomp+1
				for (j in 1:nrow(coords)) {
				compteur<-compteur+1
				if ((j == 3) & (i == bufsizes[1])) {
				esTime<-(Sys.time()-t0)*(length(bufsizes)*nrow(sites@data)/2)
				cat(paste("Time at start: ",t0,"\nApproximate run time: ", format(esTime),"\nApproximate end time: ",t0+esTime,"\n"),sep="")
				}
				if (isTRUE(all.equal(compteur/100,trunc(compteur/100)))) cat(compteur,"\n")
				flush.console()
				buff<-polycirc(i,coords[j,])
				idx<-inout(coordinates(rast),buff)
				mylist[[bufcomp]][[j]]<-rast@data[idx,att]
			}
		}

cat("Finished !\n")
return(mylist)
}
