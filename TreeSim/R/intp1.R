intp1<-function(t,lambda,mu){
	out <- 1-exp(-(lambda-mu)*t)
	out <- out /(lambda-mu*exp(-(lambda-mu)*t))
	out
	}