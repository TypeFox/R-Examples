append.df <- function(orig, add){
	app <- function(x, y){
		if(is.factor(x)){x <- as.character(x)}
		if(is.factor(y)){y <- as.character(y)}
		out <- c(x, y)
		if(is.character(out)){out <- as.factor(out)}
		out
	}
	o <- orig[,order(names(orig))]
	a <- add[,order(names(add))]
	sel <- match(names(a), names(o))
	if(sum(is.na(sel))>0){warning("Names in 'add' do not completely match")}
	sel <- sel[!is.na(sel)]
	adds <- o[,sel]
	nadds <- o[,-sel]
	out <- data.frame(matrix(NA, nrow(o)+nrow(a), ncol(adds)))
	names(out) <- names(adds)
	for(i in c(1:ncol(adds))){
		out[,i] <- app(adds[,i], a[,names(adds)[i]])
	}
	nas <- data.frame(matrix(NA, nrow(a), ncol(nadds)))
	names(nas) <- names(nadds)
	out <- cbind(out, rbind(nadds, nas))
	out <- out[,order(names(out))]
	out
}