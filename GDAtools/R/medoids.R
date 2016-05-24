medoids <- function(D,cl) {
    if(class(D)=='dist') D <- as.matrix(D)
    k <- length(table(cl))
    res <- numeric()
    for(i in 1:k) {
	x <- D[cl==i,cl==i]
	id <- 1:length(D[1,]); id <- id[which(cl==i)]
	if (length(x)>1) res[i] <- id[which.min(colMeans(x))] else res[i] <- which(cl==i)
	}
    res
    }

