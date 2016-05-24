# if vector x has NA values at the end, remove these NA elements and return the rest of the vector
rmlastna=function(x){
	while(is.na(x[length(x)])) x=x[-length(x)];
	x;}
