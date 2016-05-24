`row2col` <-
function(k){
	vec<-numeric(k*(k-1)/2)
	r<-c(0,cumsum((k-2):1))
	id<-cumsum(0:(k-1))
	for(i in 1:(k-1))
		vec[(id[i]+1):id[i+1]]<-i+r[1:i]
	vec
}

