ezglm <- function(y, x1, x2, thr = 1, family=c("gaussian","binomial")) {
	this.call <- match.call()
	family <- match.arg(family)
	thr <- as.double(thr)
	x1 <- as.double(x1)
	x2 <- as.double(x2)
	y <- as.double(y)
	no <- as.integer(length(x1))
	fit <- switch(family, 
	binomial = .Fortran("logr",no, x1, x2, y, thr, res = double(4*3)),
	gaussian = .Fortran("lsr",no, x1, x2, y, thr, res = double(4*3)))
	res = matrix(fit$res, ncol = 3)
	rownames(res) = c("(Intercept)", "x1", "x2", "x1*x2")
	colnames(res) = c("Estimate", "Std. Error", "Pr(>|t|)")
	res
} 
