probaPost <-
function(y, g, conds, pi, s, lambda) {

if(is.matrix(y) == FALSE & is.data.frame(y) == FALSE) 
	stop(paste(sQuote("y"), "must be a matrix"))
if(min(y) < 0 | sum(round(y)) != sum(y)) 
	stop(paste(sQuote("y"), "must be a matrix made up of nonnegative counts"))
if(min(rowSums(y)) == 0)
	stop(paste("at least one observation in", sQuote("y"), "contains all 0's and must be removed from the data"))
if(length(g) != 1 | g < 0 | round(g) != g) 
	stop(paste(sQuote("g"), "(the number of clusters) must be a nonnegative integer"))
if(is.vector(conds) == FALSE | length(conds) != ncol(y))
	stop(paste(sQuote("conds"), "must be a vector the same length as the number of columns in", sQuote("y")))
if(is.vector(pi) == FALSE | length(pi) != g) 
	stop(paste(sQuote("pi"), "must be a vector of length", sQuote("g")))
if(length(s) != length(conds))
	stop(paste(sQuote("s"), "and", sQuote("conds"), 
	"must be vectors the same length as the number of columns in", sQuote("y")))
if(is.matrix(lambda) == FALSE | ncol(lambda) != g | nrow(lambda) != length(unique(conds)))
	stop(paste(sQuote("lambda"), "must be a (d x g) matrix"))


n <- dim(y)[1]; cols <- dim(y)[2];
t <- matrix(0, nrow = n, ncol = g)
mean <- PoisMixMean(y, g, conds, s, lambda)
t <- matrix(unlist(lapply(1:g, function(x) .myprobafxn(k=x, y=y, pi=pi, mean=mean)), use.names=F), nrow=n, ncol=g)
## Fix problematic values of t (= 0 for all clusters)
for(j in which(rowSums(t) == 0)) {
	mean.list <- matrix(unlist(lapply(mean, function(x) x[j,]), use.names=F), ncol=length(conds), byrow=T)
	distance <- as.matrix(dist(rbind(y[j,], mean.list)))[,1]
	distance <- distance[-1]
	## If distances are exactly the same, arbitrarily pick the first as closest
	cluster <- which(distance == min(distance))[1]
	t[j,cluster] <- 1
}
## Normalize t: I think this is an error here
##t <- apply(t, 2, function(x) x/rowSums(t))
t <- t / rowSums(t)
## Smoothing prior to M-Step (0's set to 1e-10, 1's set to 1-1e-10)
epsilon <- 1e-10;maxcut <- 1-epsilon; mincut <- epsilon
t <- apply(t, 2, pmax, mincut); t <- apply(t, 2, pmin, maxcut);
## ADDED
t <- t / rowSums(t)

return(t)

}

