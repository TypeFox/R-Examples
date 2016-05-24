export <-
function(x, digits = 4, ...){
	out <- x$flux.table[,names(x$flux.table)!="all"]
	sel <- vector("logical", ncol(out))
	for(i in c(1:ncol(out))){ 
		sel[i] <- is.numeric(out[,i])
	}
	out[,sel] <- apply(out[,sel], 2, round, digits)
	write.table(out, sep="\t", row.names=FALSE, quote=FALSE, ...)
}