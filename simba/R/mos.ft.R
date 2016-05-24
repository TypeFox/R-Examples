"mos.ft" <- 
## compare a focal sampling unit (e.g. a vegetation releve)
## with the pooled surrounding units regarding the similarity 
## according to its objects (e.g. plants)
## using binary or quantitative similarity coefficients
function(x, foc=NULL, method="soerensen", quant=FALSE, binary=TRUE, ...){
	## check focal specification
	if(is.null(foc)){stop("focal plot has to be specified")}
	if(is.character(foc)){
		foc <- which(rownames(x)==foc)
	}
	sel <- x[foc,]
	x <- x[-foc,]
	if(quant){
		if(!binary){
			x.r <- colSums(x)
			x.r <- x.r/sum(x.r)
			x.c <- rbind(sel, x.r)
		}
		else{
			sel <- sel*nrow(x)
			x.c <- rbind(sel, colSums(x))
		}
		res <- vegdist(x.c, method=method)
	}
	else{
		x.c <- rbind(sel, colSums(x)>0)
		res <- sim(x.c, method=method, ...)
	}
	return(res)
}