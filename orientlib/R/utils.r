rotation.angle <- function(x) {
    x <- rotmatrix(x)@x
    apply(x,3,function(A) acos((sum(diag(A)) - 1)/2))
}

rotation.distance <- function(x, y)
{
    	x <- rotmatrix(x)@x
    	y <- rotmatrix(y)@x
    	sapply(1:dim(x)[3],
	    function(i) { 
	   	trace <- sum(x[,,i] * y[,,i])
		acos((trace - 1)/2)
       	    })
}

nearest.orthog <- function(x)
{
    d <- dim(x)
    if (length(d) < 3) d <- c(d,1)
    a <- array(x, d)
    array(apply(a,3,function(amat) {
	    sva <- La.svd(amat)
	    sva$u %*% sva$vt
       	 }),dim(a))
}

nearest.SO3 <- function(x)
{
    d <- dim(x)
    if (length(d) < 3) d <- c(d,1)
    a <- array(x, d)
    rotmatrix(array(apply(a,3,function(amat) {
	    sva <- La.svd(amat)
	    signs <- rep(1,length(sva$d))
	    if (det(amat) < 0) signs[which.min(sva$d)] <- -1
	    sva$u %*% diag(signs) %*% sva$vt
        }),dim(a)))
}
    
