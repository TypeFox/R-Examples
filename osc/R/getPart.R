cca.row <- function(pop, row){
	stopifnot(is.numeric(pop))
	stopifnot(is.matrix(pop))
	stopifnot(is.numeric(row))
	stopifnot(row<=NROW(pop))
	stopifnot(row>=0)
	the.pop <- as.integer(t(pop))
	clu <- rep(999, ncol(pop))
	out <- .C("getrow",  col=as.integer(row-1), xmax=nrow(pop), ymax=ncol(pop),  pop=as.integer(the.pop), ret=as.integer(clu))

	return(out$ret)
}
cca.col <- function(pop, col){
	stopifnot(is.numeric(pop))
	stopifnot(is.matrix(pop))
	stopifnot(is.numeric(col))
	stopifnot(col<=NCOL(pop))
	stopifnot(col>=1)
	the.pop <- as.integer(t(pop))
	clu <- rep(999, nrow(pop))
	out <- .C("getcol",  col=as.integer(col-1), xmax=nrow(pop), ymax=ncol(pop),  pop=as.integer(the.pop), ret=as.integer(clu))

	return(out$ret)
}
#cca.block <- function(pop, row,col, drow,dcol){
#	stopifnot(is.numeric(pop))
#	stopifnot(is.matrix(pop))
#	stopifnot(is.numeric(row))
#	stopifnot(row+drow<=NROW(pop))
#	stopifnot(col+dcol<=NROW(pop))
#	stopifnot(row>=1)
#	stopifnot(col>=1)
#	the.pop <- as.integer(t(pop))
#	clu <- rep(999, dcol*drow)
#	out <- .C("getblock",  row=as.integer(row-1), col=as.integer(col-1),drow=as.integer(drow),dcol=as.integer(dcol),xmax=nrow(pop), ymax=ncol(pop),  pop=the.pop, ret=as.integer(clu))
#
#	return(matrix(out$ret, ncol=dcol))
#}
