`count` <-
function(sv="species vector",x="species sample vector"){
	out=sv*0
	for (i in 1:length(sv)){
		out[i]=length((1:length(x))[x==sv[i]])
		
		}
	return(out)
	}

