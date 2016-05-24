GroupsAsText <-
function(PosCl){
	res<-array(dim=c(PosCl$dimen[1],length(PosCl$VecBin)))
	for (i in 1:PosCl$dimen[1]){
		placecounter<-1
		endcounter<-cumsum(PosCl$VecBin)
		for (j in 1:length(PosCl$VecBin)){
			res[i,j]<-paste(PosCl$group[i,placecounter:endcounter[j]],collapse="")
			placecounter<-placecounter+PosCl$VecBin[j]
			}
		}
	return(res)	
	}

