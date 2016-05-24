"sticks" <-
function (num=100, color="black", lwd=1, seed=NULL) 
{
	if(length(seed)) set.seed(seed)

	ends <- array(NA, c(num, 4))
	onetwo <- 1:2
	threefour <- 3:4
	ends[1,onetwo] <- runif(2)
	ends[1,threefour] <- ends[1, onetwo] + runif(2, -1, 1)
	for(i in 2:num) {
		ends[i, onetwo] <- .5 * (ends[i-1, onetwo] + 
			ends[i-1, threefour])
		ends[i, threefour] <- ends[i, onetwo] + runif(2, -1, 1)
	}
	plot(.5, .5, xlim=range(ends[, c(1,3)]), ylim=range(ends[, c(2,4)]),
		type="n", xlab="", ylab="", axes=FALSE)
	segments(ends[,1], ends[,2], ends[,3], ends[,4],
	         col=safesample(color), lwd=safesample(lwd))
}

