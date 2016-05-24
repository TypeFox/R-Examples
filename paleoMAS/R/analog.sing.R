analog.sing <-
function(fossil,base,age,dca=FALSE,method="euclidean",
	binary=FALSE)
{
library(vegan)
{
	distance<-matrix(nrow=nrow(fossil),ncol=2)
	colnames(distance)<-c("Age","Distance")
	distance[,1]<-age
	if(dca==TRUE){
		decorana(fossil)->dca.f
		scores(dca.f,display="sites")->sscores
		as.matrix(vegdist(sscores,method=method,
			binary=binary))[,base]->distance[,2]
		}
	else{
		if(method=="schord"){
			for(i in 1:nrow(fossil)){
				distance[i,2]<-sum(((fossil[i,
					]^0.5)-(fossil[base,]^0.5))^2)
				}
		}
		else{
			as.matrix(vegdist(fossil,method=method,
				binary=binary))[,base]->distance[,2]			}
		}
	plot(distance[,c(2,1)],ylim=c(max(age),min(age)),
			type="l")
}
return(distance)
}

