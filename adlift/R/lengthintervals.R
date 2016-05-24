`lengthintervals` <-
function (X, I, type = "midpoints", neighbours, closest) 
{

o <- as.column(order(X))
lengths <- matrix(0, 1, length(X))
d <- neighbours

if (closest) {
	initialnbrs <- matrix(0, length(X), neighbours)
        if (type == "average") {
        	for (i in 1:length(X)) {
                	out <- getnbrs(X, o[i], order(X), neighbours, closest)
                	initialnbrs[i, ] <- out$nbrs
                	lengths[i] <- sum(abs(rep(X[o[i]], times = d) - X[initialnbrs[i, ]]))/d
            	}
        }
	if (type == "midpoints") {
		lengths<-I[2:length(I)]-I[1:length(X)] 
        }
}
else {
	lengths<-I[2:length(I)]-I[1:length(X)]
}
    
return(as.row(lengths))
}

