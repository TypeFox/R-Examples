rsse<-function(a,b,v=1){

if(!is.matrix(a)){
	a<-as.matrix(a)
}

b<-matrix(b,byrow=T,nrow=nrow(a),ncol=ncol(a))
v<-matrix(v,byrow=T,nrow=nrow(a),ncol=ncol(a))

stda<-((a-b)^2)/v 

rsse<-sqrt(sum(rowMeans(stda)))

rsse


}
