maxstat.table<-function(x,...) {
 if (dim(x)[1]> 2) x<-t(x)
 if (dim(x)[1]> 2) stop("table should be 2x3")
 if (dim(x)[2]==1) warning("this SNP is monomorphic")
 maxstat.matrix(x, ...)
 }
