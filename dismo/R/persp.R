
if (!isGeneric("persp")) {
	setGeneric("persp", function(x,...)
		standardGeneric("persp"))
}	


setMethod("persp", signature(x='DistModel'), 
function(x, var1, var2, at=median, ...) {
	
	d <- x@presence
	
	if (missing(var1)) { var1 <- colnames(d)[1] }
	if (missing(var2)) { var2 <- colnames(d)[2] }
	xi <- which(colnames(d)==var1)
	yi <- which(colnames(d)==var2)
	if (length(xi) == 0) { stop('x not found') }
	if (length(yi) == 0) { stop('y not found') }
	
	v1 <- d[, var1]
	v2 <- d[, var2]
	if (is.factor(v1) | is.factor(v2)) {
		stop('I cannot do this for factors yet')
	} 

	v1 <- range(v1)
	v2 <- range(v2)
	v1 <- v1[1] + -5:55 * 2 * (v1[2]-v1[1])/100
	v2 <- v2[1] + -5:55 * 2 * (v2[2]-v2[1])/100
	v <- cbind(rep(v1, each=61), rep(v2, 61))
		
	m <- as.numeric(apply(d[,-c(xi,yi)], 2, at))
	m <- data.frame(matrix(rep(m, each=nrow(v)), nrow=nrow(v)))
	colnames(m) <- colnames(d)[-c(xi, yi)]
	a <- cbind(v, m)
	colnames(a)[1:2] <- c(var1,var2)
	p <- predict(x, a)

	z <- matrix(p, ncol=61)
	persp(x=v1, y=v2, z=z, ...)
}
)

	
