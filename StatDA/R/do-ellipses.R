do.ellipses <- function (acov, pos, ...)		## internal function
{
# Plot ellipses
	acov = cov2cor (acov)
	cov.svd <- svd(acov, nv = 0)
	r <- cov.svd[["u"]] %*% diag(sqrt(cov.svd[["d"]]))

	m <- 100

	alphamd <- c(1/3)

    e1md <- cos(c(0:m)/m * 2 * pi) * alphamd
    e2md <- sin(c(0:m)/m * 2 * pi) * alphamd
    emd <- cbind(e1md, e2md)
    ttmd <- t(r %*% t(emd)) + rep(1, m + 1)
    lines(ttmd[, 1] + pos[1], ttmd[, 2]+ pos[2], ...)
}
