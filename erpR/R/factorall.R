factorall <-
function(x){
	if (length(dim(x))>1){
		for(i in 1:length(x)){
			{
				x[,i]=factor(x[,i])}
			}
		return(x)
		} else {
		x=factor(x)
		return(x)
		}
	}
