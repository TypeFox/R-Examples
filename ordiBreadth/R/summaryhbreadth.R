summaryhbreadth <-
function(specialization,round=5,do.order=FALSE,by="Richness"){
	res<-data.frame(cbind(specialization$species,apply(specialization$group.vectors==TRUE,1,sum),round(specialization$tot.breadth,digits=round),round(specialization$scaled.breadth,digits=round)))
	colnames(res)<-c("Herbivore","Richness","Breadth","ScaledBreadth")
	if(do.order==TRUE){
		neworder<-order(as.numeric(as.vector(res[,colnames(res)==by])))
		res<-res[neworder,]
		}	
	return(res)
	}
