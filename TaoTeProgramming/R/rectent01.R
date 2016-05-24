"rectent01" <-
function (x=array(sample(c(TRUE, FALSE), 4e6, replace=TRUE), c(2000,2000)),
	iterations=100, inner=20, verbose=TRUE) 
{
	initEnt <- entropy(x)
	if(verbose) cat(initEnt, "full matrix", date(), "\n")
	dx <- dim(x)
	jseq <- 1:inner
	choice <- c(TRUE, FALSE)
	for(i in 1:iterations) {
		xout <- sort(sample.int(dx[1], 2))
		yout <- sort(sample.int(dx[2], 2))
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

