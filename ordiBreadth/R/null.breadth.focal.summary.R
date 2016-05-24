null.breadth.focal.summary <-
function(null.breadth.focal.out,quantiles=c(0.025,0.975),round=5,scaled=FALSE){
	species<-ODB<-scale.factor<-richness<-adj.richness<-lq<-uq<-NA
	for(i in 1:length(null.breadth.focal.out)){
		species[i]<-null.breadth.focal.out[[i]]$species
		ODB[i]<-null.breadth.focal.out[[i]]$observed.breadth
		scale.factor[i]<-null.breadth.focal.out[[i]]$scale.factor
		richness[i]<-null.breadth.focal.out[[i]]$totalplantrichness
		adj.richness[i]<-null.breadth.focal.out[[i]]$modplantrichness
		lq[i]<-quantile(null.breadth.focal.out[[i]]$null,quantiles[1])
		uq[i]<-quantile(null.breadth.focal.out[[i]]$null,quantiles[2])
			}
			if(scaled==TRUE){
				ODB<-round(ODB/scale.factor,round);lq<-round(lq/scale.factor,round);uq<-round(uq/scale.factor,round)
			}else{ODB<-round(ODB,round);lq<-round(lq,round);uq<-round(uq,round)}
			scale.factor<-round(scale.factor,round)
		res<-data.frame(cbind(species,richness,adj.richness,ODB,lq,uq,scale.factor))
		return(res)
	}
