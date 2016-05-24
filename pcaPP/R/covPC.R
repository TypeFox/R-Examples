covPC <- function (x, k = ncol (x$loadings), method)
{
	if (!any(class(x) == "princomp"))
		stop ("Invalid parameter \x22k\x22: Data type princomp expected!")
	if (length (x$sdev) != length (x$center))
		warning ("Calculating a rank", length (x$sdev), "- covariance matrix")

	ret = list()
	k = min (ncol (x$loadings), k)

	ret$cov = x$loadings[,1:k] %*% diag (x$sdev[1:k]^2) %*% t(x$loadings[,1:k])
	ret$center = x$center
	if (missing (method))
		ret$method = "Covariance estimation based on PCs"
	else
		ret$method = method

	class (ret) <- "covPC"

	return (ret)
} 

covPCAgrid <- function (x, control)
{
	pcs = PCAgrid (x, k = ncol(x), control = control)

	ret = list()
	ret$cov = pcs$loadings %*% diag (pcs$sdev^2) %*% t(pcs$loadings)
	ret$center = pcs$center
	ret$method = "Robust cov - estimation based on PCs (grid mode)"

	if (!missing (control) && !is.null (control$method))
		ret$method = paste ("Robust cov - estimation based on PCs (grid mode -", control$method, ")", sep = "")
	else
		ret$method = "Robust cov - estimation based on PCs (grid mode)"

	class (ret) <- "covPC"

	return (ret)
}

covPCAproj <- function (x, control)
{
	pcs = PCAproj (x, k = ncol(x), control = control)

	ret = list()
	ret$cov = pcs$loadings %*% diag (pcs$sdev^2) %*% t(pcs$loadings)
	ret$center = pcs$center

	if (!missing (control) && !is.null (control$method))
		ret$method = paste ("Robust cov - estimation based on PCs (projection mode - ", control$method, ")", sep = "")
	else
		ret$method = "Robust cov - estimation based on PCs (projection mode)"

	class (ret) <- "covPC"

	return (ret)
}
