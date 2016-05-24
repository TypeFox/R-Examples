negent2D <-
function(y,m=100) {
	thetas <- (1:m)*pi/m
	ngnt <- NULL
	for(theta in thetas) {
		x <- y%*%c(cos(theta),sin(theta))
		ngnt <- c(ngnt,negent(x))
	}
	i <- imax(ngnt)
	g <- c(cos(thetas[i]),sin(thetas[i]))
	g <- cbind(g,c(-g[2],g[1]))
	list(vectors = g,values = c(ngnt[i],negent(y%*%g[,2])))
}
