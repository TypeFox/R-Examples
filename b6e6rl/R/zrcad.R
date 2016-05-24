zrcad <-
function(y,a,b){
zrc <- which((y<a) | (y>b))
if (is.empty(zrc) != TRUE){
for (i in 1:length(zrc)){
	while((y[i]<a[i]) | (y[i]>b[i])){
		if (y[i]>b[i])
			y[i] <- 2*b[i]-y[i]
		if(y[i]<a[i])
			y[i] <- 2*a[i]-y[i]			
	}
}
}
result <- y
return(result)
}
