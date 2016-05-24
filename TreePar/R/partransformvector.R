partransformvector <-
function(p){
	numbs <- (length(p)+1)/3
	l<-vector()
	mu<- vector()
	t<-vector()
	for (i in 1:numbs){
		l<-c(l, p[i+numbs]/(1-p[i]))
    	mu<- c(mu,p[i] * p[i+numbs]/(1-p[i]))  
		if (i<numbs){
			t<-c(t,p[i+2*numbs])
			}
		}
	t<-c(0,t)
	rbind(l,mu,t)
	}

