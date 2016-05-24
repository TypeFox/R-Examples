samplesize.d <- function(delta, power=.80, sig.level=.05){
	n <- 2*((qnorm(sig.level/2,lower.tail=FALSE)-qnorm(power,lower.tail=FALSE))/
	delta)^2+(qnorm(sig.level/2,lower.tail=FALSE)^2)/4
	return(ceiling(n))
}