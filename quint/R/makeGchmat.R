makeGchmat <-
function(gm,varx,splitpoint){
  #creates an indicator matrix of the child nodes
	#gm = indicator vector of persons in the parent node that is split
	#varx =  unsorted datavector of splitting predictor X
	#splitpoint is optimal threshold on varx used for splitting
	Gchmat<-matrix(nrow=length(gm),ncol=2)
	Gchmat<-cbind(ifelse(gm==1&varx<=splitpoint,1,0),ifelse(gm==1&varx>splitpoint,1,0))
	return(Gchmat)}
