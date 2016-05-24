summaryordi <-
function(ordi.out,round=5){
	#ordi.out is an object from ordi.focal.drop
	species<-NA
	ODB<-NA
	focal.breadth<-NA
	scaled.ODB<-NA
	scaled.focal.breadth<-NA
		for (i in 1:length(ordi.out)){
			species[i]<-ordi.out[[i]]$species
			ODB[i]<-round(ordi.out[[i]]$ODB,digits=round)
			focal.breadth[i]<-round(ordi.out[[i]]$focal.breadth,digits=round)
			scaled.ODB[i]<-round(ordi.out[[i]]$ODB.scaled,digits=round)
			scaled.focal.breadth[i]<-round(ordi.out[[i]]$focal.scaled.breadth,digits=round)
			}
		res<-data.frame(cbind(species,ODB,focal.breadth,scaled.ODB,scaled.focal.breadth))	
	return(res)
	}
