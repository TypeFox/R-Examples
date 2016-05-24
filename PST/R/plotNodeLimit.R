##
plotNodeLimit <- function (x1, x2, subtree, max.level) {
	inner <- !all(subtree@leaf, na.rm=TRUE) && x1 != x2 && !subtree@order==max.level

	if (inner) {
	        ## K <- length(subtree)
		K <- which.child(subtree)
        	mTop <- summary(subtree, max.level=max.level, segmented=FALSE)@leaves
        	limit <- integer(length(K))
		names(limit) <- K
        	xx1 <- x1

        	for (k in K) {
            		m <- summary(subtree[[k]], max.level=max.level, segmented=FALSE)@leaves
            		xx1 <- xx1 + ((x2 - x1) * m/mTop)
			limit[k] <- xx1
		}
		limit <- c(x1, limit)
	} else {
		limit <- c(x1, x2)
	}
	## mid <- attr(subtree, "midpoint")
	mid <- length(subtree)/2

    x <- mean(c(x1, x2))
    list(x = x, limit = limit)
}

