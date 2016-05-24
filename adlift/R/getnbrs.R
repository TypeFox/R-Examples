`getnbrs` <-
function (X, remove, pointsin, neighbours, closest) 
{

N <- length(pointsin)
d <- neighbours
nbrs <- distances<-NULL
r <- which(remove == pointsin)
range <- 2:(N - 1)

if (is.element(r, range)) {
	checkl <- (r - (1:d) < 1)
        leftneigh <- d - sum(checkl)
        checkr <- (r + (1:d) > N)
        rightneigh <- d - sum(checkr)
        nbrs <- pointsin[r - leftneigh + (1:leftneigh) - 1]
        nbrs[(1:rightneigh) + leftneigh] <- pointsin[r + (1:rightneigh)]
}
if (!closest) {
        if (r == 1) {
        	nbrs <- pointsin[2]
        	leftneigh <- 0
        	rightneigh <- 1
        }
        if (r == N) {
        	nbrs <- pointsin[N - 1]
        	leftneigh <- 1
        	rightneigh <- 0
        }
        index <- setdiff((r - leftneigh):(r + rightneigh), r)
}
else {
        if (r == 1) {
        	index <- r + 1
        	nbrs <- pointsin[index]
        	leftneigh <- 0
        	rightneigh <- 1
        }
        if (r == N) {
        	index <- (r - 1)
        	nbrs <- pointsin[index]
        	leftneigh <- 1
        	rightneigh <- 0
        }
        if (is.element(r, range)) {
        	distances <- abs(X[remove] - X[pointsin[r - leftneigh + (1:leftneigh) - 1]])
        	distances[(1:rightneigh) + leftneigh] <- abs(X[remove] - X[pointsin[r + (1:rightneigh)]])
        	d1 <- min(d, N - 1)
        	q <- order(distances)[1:d1]
        	nbrs <- nbrs[q]
        	index <- setdiff((r - leftneigh):(r + rightneigh), r)
        	index <- index[q]
        	B <- (index[q[1:d1]] < r)
        	numleft <- sum(B)
        	leftneigh <- numleft
        	rightneigh <- d1 - leftneigh
        }
}

index <- sort(index)
nbrs <- as.vector(pointsin[index])

return(list(nbrs = nbrs, index = index))
}

