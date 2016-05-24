B2M <-
function(x) {
	x[x==0] <- min(x[x!=0])
	x[x==1] <- max(x[x!=1])
	log(x / (1-x))
}
