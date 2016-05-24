Rmatsolve<-function(m){

if(!is.matrix(m)){
	m<-as.matrix(m)
}

if (nrow(m)==1){     #since diag doesn't like 1x1 matrices
	inv<-1/m
}
else{
	e <- eigen(m, symmetric=TRUE)
	ev <- e$vectors
	inv<-ev %*% diag(1/e$values) %*% t(ev)
}

inv
}
