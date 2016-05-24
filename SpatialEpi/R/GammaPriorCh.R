GammaPriorCh <-
function(theta, prob, d){
	a <- d/2
	b <- 0.5*2*a*theta^2/qt(p=prob,df=2*a)^2
	cat("Gamma Parameters: ",a,b,"\n")
	list(a=a,b=b)
}
