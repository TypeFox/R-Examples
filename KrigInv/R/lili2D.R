lili2D <- function(x){
	x1 <- (x[,1])*14-7
	x2 <- (x[,2])*14-7
	f1 <- 3 + 0.1*(x1-x2)^2 - (x1+x2)/sqrt(2)
	f2 <- 3 + 0.1*(x1-x2)^2 + (x1+x2)/sqrt(2)
	f3 <- (x1 - x2) + 6/sqrt(2)
	f4 <- -(x1 - x2) + 6/sqrt(2)
	f<- pmin(f1,f2,f3,f4)
	return(f)
}
