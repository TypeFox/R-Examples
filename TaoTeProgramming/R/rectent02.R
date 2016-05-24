"rectent02" <-
function (x=array(sample(c(TRUE, FALSE), 4e6, replace=TRUE), c(2000,2000)),
	iterations=100, inner=20, fraction=c(.01, .1), seed=NULL, verbose=TRUE) 
{
	if(length(seed)) set.seed(seed)
	initEnt <- entropy(x)
	if(verbose) cat(initEnt, "full matrix", date(), "\n")
	dx <- dim(x)
	jseq <- 1:inner
	choice <- c(TRUE, FALSE)
	xrange <- round(dx[1] * min(fraction)): round(dx[1] * max(fraction))
	yrange <- round(dx[2] * min(fraction)): round(dx[2] * max(fraction))
	xtop <- dx[1] - max(xrange)
	ytop <- dx[2] - max(yrange)
	for(i in 1:iterations) {
		xout <- sample.int(xtop, 1) + c(0, sample(xrange, 1))
		yout <- sample.int(ytop, 1) + c(0, sample(yrange, 1))
		ixmat <- x[ xout[1]:xout[2], yout[1]:yout[2] ]
		innerEnt <- entropy(ixmat)
		if(verbose) cat(innerEnt, "initial iteration", i, date(), "\n")
		ixlen <- length(ixmat)
		for(j in jseq) {
			newx <- xor(ixmat, sample(choice, ixlen, replace=TRUE))
			newEnt <- entropy(newx)
			if(newEnt < innerEnt) {
				innerEnt <- newEnt
				ixmat <- newx
				if(verbose) {
					cat(" ", innerEnt, "at inner", j, "\n")
				}
			}
		}
		x[ xout[1]:xout[2], yout[1]:yout[2] ] <- ixmat
	}
	finEnt <- entropy(x)
	if(verbose) cat(finEnt, "full matrix", date(), "\n")
	x
}

