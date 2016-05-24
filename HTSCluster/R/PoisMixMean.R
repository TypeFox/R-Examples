PoisMixMean <-
function(y, g, conds, s, lambda) {
if(is.matrix(y) == FALSE & is.data.frame(y) == FALSE) 
	stop(paste(sQuote("y"), "must be a matrix"))
if(min(y) < 0 | sum(round(y)) != sum(y)) 
	stop(paste(sQuote("y"), "must be a matrix made up of nonnegative counts"))
if(min(rowSums(y)) == 0)
	stop(paste("at least one observation in", sQuote("y"), "contains all 0's and must be removed from the data"))
if(length(g) != 1)
	stop(paste(sQuote("g"), "(the number of clusters) must be a nonnegative integer"))
if(g < 0 | round(g) != g) 
	stop(paste(sQuote("g"), "(the number of clusters) must be a nonnegative integer"))
if(is.vector(conds) == FALSE | length(conds) != ncol(y))
	stop(paste(sQuote("conds"), "must be a vector the same length as the number of columns in", sQuote("y")))
if(length(s) != length(conds))
	stop(paste(sQuote("s"), "and", sQuote("conds"), 
	"must be vectors the same length as the number of columns in", sQuote("y")))
if(is.matrix(lambda) == FALSE | ncol(lambda) != g | nrow(lambda) != length(unique(conds)))
	stop(paste(sQuote("lambda"), "must be a (d x g) matrix"))
n <- dim(y)[1];cols <- dim(y)[2];
r <- as.vector(table(conds))
w <- rowSums(y)
mean.mat <- vector("list", g)
w.mat <- matrix(rep(w, times = cols), nrow = n, ncol = cols)
s.mat <- matrix(rep(s, each = n), nrow = n, ncol = cols)
mean.mat <- lapply(1:g, function(x)
	.myloopfxn(x, lambda=lambda, w.mat=w.mat, s.mat=s.mat, r=r, n=n, cols=cols))
return(mean.mat)
}

