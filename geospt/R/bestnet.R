bestnet <- function(optimnet){
	if(class(optimnet)!= "rbga") stop(paste("object must be of class rbga"))
	j <- which.min(optimnet$evaluations)
	a <- optimnet$population[j,]
	n <- dim(optimnet$population)[2]/2
	b <- cbind(a[1:n],a[(n+1):(2*n)])
	p.c <- as.data.frame(round(b, 1))
	names(p.c) <- c("x", "y")
	coordinates(p.c) = c("x", "y")
	p.c
}
