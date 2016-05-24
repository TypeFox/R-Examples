index.mat.generation <-
function(word.vec,pos.vec){

	mat.index<-matrix(0,nrow=length(word.vec),ncol=length(pos.vec),dimnames=list(word.vec,pos.vec))
	flat.index<-cbind(1:length(mat.index),rep(word.vec,length(pos.vec)),rep(pos.vec,each=length(word.vec)))

	return(list(mat=mat.index,flat=flat.index))
}
