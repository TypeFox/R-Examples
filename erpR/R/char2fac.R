char2fac <-
function(x){
		for (i in 1:length(x)) {
			if (is.character(x[,i])) {
					x[,i]=factor(x[,i])
			}
		}
return(x)
}
