`compVarLevs` <-
function(listX){
	listX<-lapply(listX,function(x) rowSums(x)>0)
	catX<-matrix(unlist(listX),ncol=length(listX))
	rowSums(catX)
}

