
which.suck=function(vectors, dim){
 if(dim>1){which(vectors<0.10 | vectors>0.4)
	}else{which(vectors<0.30 | vectors>0.55)}

}
