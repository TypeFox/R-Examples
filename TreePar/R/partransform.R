partransform <-
function(a,r){
	l<- r/(1-a)
	mu <- a*l
	c(l,mu)
	}

