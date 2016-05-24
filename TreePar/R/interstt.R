interstt <- function(time,t) {
	res<-1
	j=0
	i=1
	while (j==0) {
		if (i==length(t)){
			j=1
		} else {
			if (time<=t[i+1]) {j=1} else {i<-i+1}
		}
	}
	i
	}	
