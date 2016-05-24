`removeRedundancy` <-
function(mat,vec.ia,red=0){
	mat<-as.data.frame(mat)
	for(i in 1:length(vec.ia))
		vec.ia[i]<-evalRedundancy(mat,vec.ia[i],red=red)
	vec.ia
}

