`dist2Mats` <-
function(mat1,mat2,method,p=2){
	n.row<-nrow(mat1)
	mat<-rbind(mat1,mat2)
	vec.dist<-dist(mat,method=method,p=p)
	mat<-as.matrix(vec.dist)
	mat[1:n.row,(n.row+1):nrow(mat)]
}

