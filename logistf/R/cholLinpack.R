cholLinpack <-
function(x=diag(3)) {
if(nrow(x) != ncol(x))
stop("wrong dimensions!")
res <- .C("linpack_choleski", res=as.double(x), 
as.integer(nrow(x)), PACKAGE="logistf")$res
matrix(res, nrow(x))
}

