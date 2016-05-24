"rm.outlier" <-
function (x, fill = FALSE, median = FALSE, opposite = FALSE) 
{
    if (is.matrix(x)) 
        apply(x, 2, rm.outlier, fill = fill, median = median, opposite = opposite)
    else if (is.data.frame(x)) 
        as.data.frame(sapply(x, rm.outlier, fill = fill, median = median, opposite = opposite))
    else {
	res <- x
	if (!fill) res[-which(x == outlier(x,opposite))]
	else {
		if (median) res[which(x == outlier(x,opposite))]<-median(x[-which(x == outlier(x,opposite))])
		else res[which(x == outlier(x,opposite))]<-mean(x[-which(x == outlier(x,opposite))])
	res
	}
}
}

