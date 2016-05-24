best.band <- function(sam, bandwidth, cv.fold, norm){
	p <- dim(sam)[2]
	n <- dim(sam)[1]
	fold.size <- round(n/cv.fold)
	sam.idx <- sample(1:n, size = n, replace = FALSE)
	n.bandwidth <- length(bandwidth)
	diff.norm <- matrix(0, cv.fold, n.bandwidth)
	for(i in 1:cv.fold){
		if(i == cv.fold){
			temp.idx <- sam.idx[((i - 1)*fold.size + 1):n]
		}else{
			temp.idx <- sam.idx[((i - 1)*fold.size + 1):(i*fold.size)]
		}
		sam.train <- sam[-temp.idx,]
		sam.test <- sam[temp.idx,]
		for(j in 1:n.bandwidth){
			sam.train.cov <- cov(sam.train)
			sam.train.cov.band <- sam.train.cov
			sam.train.cov.band[abs(row(sam.train.cov.band) - col(sam.train.cov.band)) > bandwidth[j]] <- 0
			sam.test.cov <- cov(sam.test)
			diff.norm[i, j] <- norm(sam.train.cov.band - sam.test.cov, type = norm)
		}
	}
	diff.norm <- colMeans(diff.norm)
	best <- which(diff.norm == min(diff.norm))
	best <- best[1]
	return(bandwidth[best])
}