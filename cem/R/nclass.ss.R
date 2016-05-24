#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2006 Author Hideaki Shimazaki
# Department of Physics, Kyoto University
# shimazaki at ton.scphys.kyoto-u.ac.jp
# Please feel free to use/modify/distribute this program.
nclass.ss <- function(x){
	N <- 2:100
	C <- numeric(length(N))
	D <- C
	for (i in 1:length(N)) {
		D[i] <- diff(range(x))/N[i]
		edges = seq(min(x),max(x),length=N[i])
		hp <- hist(x, breaks = edges, plot=FALSE )
		ki <- hp$counts
		k <- mean(ki)
		v <- sum((ki-k)^2)/N[i]
		C[i] <- (2*k-v)/D[i]^2	#Cost Function
	}
	idx <- which.min(C)
    N[idx]
#	optD <- D[idx]
#
#	edges <- seq(min(x),max(x),length=N[idx])
#	h = hist(x, breaks = edges )
#	rug(x)
#
#	return(h)
#    edges
}