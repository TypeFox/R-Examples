IndexModify <-
function(index,obj.index,...){
	c <- ncol(index)
	id <- NULL
	for(i in 1:c){
		for(j in 1:c){
			k.id <- which(index[,i] == obj.index[j])
			id <- union(id,k.id)
		}
	}
	if(length(id)>0)index <- matrix(index[-id,],,c)
	return(index)
}
