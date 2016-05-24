f.post.diff <- function(coeff, covar){
##
## FOR THE HAPLOTYPE FREQUENCIES, SUBTRACT FIRST PARAMETER FROM THE REST,
## TO "NORMALIZE". DO THIS BY MULTIPLYING WITH DIFFERENCE MATRIX
.names <- names(coeff[[1]])
#
.D <- diag(length(coeff[[1]]))
dimnames(.D) <- list(.names, .names)
#
.tmp <- f.coefnames(.names)
.mf <- .tmp$haplo.freq
.mf1 <- match(.mf[1], .names) # not necessarily certain that mf is at the start
.D[.mf, .mf[1]] <- -1
.D <- .D[-.mf1,]
#
.coef <- lapply(coeff, function(x) .D %*% x)
.cov <- lapply(covar, function(x) .D %*% x %*% t(.D))
.vis <- F
if(.vis){
	f.vis(.D, vis = .vis)
	f.vis(coeff[[1]], vis = .vis)
	f.vis(.D %*% coeff[[1]], vis = .vis)
	#
	f.vis(covar, vis = .vis)
	f.vis(.coef, vis = .vis)
	f.vis(.cov, vis = .vis)
}
#
##
return(list(coeff = .coef, covar = .cov))
}
