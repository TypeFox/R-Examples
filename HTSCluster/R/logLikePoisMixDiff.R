logLikePoisMixDiff <-
function(y, mean.new, pi.new, mean.old, pi.old) {
if(is.matrix(y) == FALSE & is.data.frame(y) == FALSE) 
	stop(paste(sQuote("y"), "must be a matrix"))
if(min(y) < 0 | sum(round(y)) != sum(y)) 
	stop(paste(sQuote("y"), "must be a matrix made up of nonnegative counts"))
if(min(rowSums(y)) == 0)
	stop(paste("at least one observation in", sQuote("y"), "contains all 0's and must be removed from the data"))
if(length(mean.old) != length(pi.old)) 
	stop(paste(sQuote("mean.old"), "must be a list of the same length as", sQuote("pi.old")))
if(length(mean.new) != length(pi.new)) 
	stop(paste(sQuote("mean.new"), "must be a list of the same length as", sQuote("pi.new")))
if(length(mean.old) != length(mean.new)) 
	stop(paste(sQuote("mean.new"), "must be a list of the same length as", sQuote("mean.old")))

n <- dim(y)[1]; cols <- dim(y)[2]
g <- length(pi.old)
y <- matrix(y, nrow = n, ncol = cols)
num.tmp <- matrix(unlist(lapply(mean.new, function(x) .myfxn(var1=y, var2=x)), use.names=F, recursive=F), nrow=n, ncol=g)
den.tmp <- matrix(unlist(lapply(mean.old, function(x) .myfxn(var1=y, var2=x)), use.names=F, recursive=F), nrow=n, ncol=g)
num <- rowSums(matrix(pi.new * matrix(exp(-num.tmp),byrow=T,nrow=g),byrow=T,ncol=g,nrow=n))
den <- rowSums(matrix(pi.old * matrix(exp(-den.tmp),byrow=T,nrow=g),byrow=T,ncol=g,nrow=n))
num <- ifelse(num == 0, NA, num)
den <- ifelse(den == 0, NA, den)
ll.tmp <- ifelse(is.nan(log(num) - log(den)) == TRUE, NA, log(num) - log(den))
ll <- sum(ll.tmp, na.rm = TRUE)
## In case one or both of the log-likelihoods is infty or -infy, set difference to large number
if(ll == -Inf | ll == Inf | is.nan(ll) == TRUE) {ll = 100}
return(ll)


}

