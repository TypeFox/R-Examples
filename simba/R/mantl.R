"mantl" <-
function(x, y, method = "pearson", permutations = 1000, sub = NULL, loop = FALSE, ...){
	p <- permutations
	x <- as.matrix(x)
	y <- as.matrix(y)
	sel <- row(x)>col(x)
	if(!is.null(sub)){ sel <- (sel+sub)==2 }
	stat <- cor(x[sel], y[sel], method = method, ...)
	if(loop){
		when <- seq(0, permutations, floor(permutations/10))
		perms <- rep(NA, permutations)
		for(i in c(1:permutations)){
			py <- sample(c(1:nrow(x)))
			perms[i] <- cor(x[sel], y[py,py][sel], method = method, ...)
			if(sum(when==i)==1) cat(".",sep="")
		}
	}
	else{
		permsy <- lapply(1:p, function(z) sample(c(1:nrow(x))))
		perms <- sapply(permsy, function(z) cor(x[sel], y[z,z][sel], method = method, ...))
	}
	if (stat >= 0) {
		signif <- sum(perms >= stat)/p
	}
	else {
		signif <- sum(perms <= stat)/p
	}
	N <- length(x[sel])
	res <- list(call=match.call(), method=method, statistic=stat, signif=signif, n=N, permutations=p, perms=perms)
	class(res) <-"permcor"
	res
}