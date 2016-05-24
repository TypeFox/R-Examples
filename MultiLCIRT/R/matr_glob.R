matr_glob <-
function(l,type="g"){

# derive matrices
	Co = cbind(-diag(l-1),diag(l-1))
	if(type=="g"){
		Ma = cbind(lower.tri(matrix(1,l-1,l-1),diag=T),rep(0,l-1));
		Ma = rbind(Ma,1-Ma)
	}
	if(type=="l"){
		Ma = rbind(cbind(diag(l-1),rep(0,l-1)),
        		       cbind(rep(0,l-1),diag(l-1)))
	}
	out = list(Co=Co,Ma=Ma)
}
