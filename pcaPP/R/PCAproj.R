
PCAproj <- function (x, k = 2, method = c ("mad", "sd", "qn"), CalcMethod = c("eachobs", "lincomb", "sphere"), nmax = 1000, update = TRUE, scores = TRUE, maxit = 5, maxhalf = 5, scale = NULL, center = l1median_NLM, zero.tol = 1e-16, control)
{
	if (!missing (control))
		###ParseControlStructure (control, c("k", "method", "CalcMethod", "nmax", "update", "scores", "maxit", "maxhalf"))
		.ParseControlStructure (control, c("k", "method", "CalcMethod", "nmax", "update", "scores", "maxit", "maxhalf", "center", "scale", "zero.tol"))

	method <- .getScaleMethod (method)

	CalcMethod <- match.arg (CalcMethod[1], c("eachobs", "lincomb", "sphere"))

	x <- .Conv2Matrix (x, substitute (x))

	n = nrow (x)
	p = ncol (x)

	if( k > min(n,p))
		stop ('k too large')

	if(p > n)
	{
			svdx = svd(t(x))
			x = svdx$v %*% diag(svdx$d)
			pold=p
			p=n
	} else
			pold=p

	DataObj = ScaleAdv (x, scale = scale, center = center)

    if (pold > n) # center and scale must have original data dimension:
	{
		DataObj$center <- as.vector(svdx$u%*%DataObj$center)
		DataObj$scale <- ScaleAdv(x%*%t(svdx$u),center=NULL,scale=scale)$scale
	}

	y = DataObj$x

#	m = l1median(x)
#	y = t(t(x) - m)

	if (scores)
		scoresize <- n * k
	else
		scoresize <- 0

	if (CalcMethod == "lincomb")
	{
		update = FALSE
		if (nmax > n)
		{
			aux = matrix (runif ((nmax-n) * n), nrow = nmax-n)
			##y = rbind (y, t(t(aux %*% x) - DataObj$center))
			y <- rbind (y, aux %*% y) ## use this instead?
		}
	}
	else if (CalcMethod == "sphere")
	{
		update = FALSE
		if(nmax  >n)
			#y[(n+1):nmax,] =	rmvnorm(nmax-n, rep(0,p), diag (p))
			y = rbind (y,		rmvnorm(nmax-n, rep(0,p), diag (p)))
	}

	nn = nrow (y)

	if (update)
		ret.C = .C ("pcaProj_up", PACKAGE="pcaPP"
				, as.integer (c(nn, p, n, k, method, scores, maxit, maxhalf))
				, as.double (zero.tol)
				, as.double (y)
				, scores = double (scoresize)
				, loadings = double (p * k)
				, lambda = double (k)
			)
	else
		ret.C = .C ("pcaProj",  PACKAGE="pcaPP"
				, as.integer (c(nn, p, n, k, method, scores))
				, as.double (zero.tol)
				, as.double (y)
				, scores = double (scoresize)
				, loadings = double (p * k)
				, lambda = double (k)
			)

	veig = matrix (ret.C$loadings, ncol = k)

	idx.mo <- ret.C$lambda == -1
	if (any (idx.mo))
	{
		veig [, idx.mo] <- .Null (veig[, !idx.mo])
		ret.C$lambda[idx.mo] <- 0
	}

   if(pold>n)
		veig = svdx$u %*% veig

	if (scores)
		.DataPostProc (DataObj, ret.C$lambda, veig, matrix (ret.C$scores, ncol = k) , match.call(), scores)
	else
		.DataPostProc (DataObj, ret.C$lambda, veig, NULL, match.call(), scores)
}

.Null <- function (M)	##	 Null function from package MASS -> 2do: move calls to LAPACK to C++ code
{
    tmp <- qr(M)
    set <- if (tmp$rank == 0L)
        1L:ncol(M)
    else -(1L:tmp$rank)
    qr.Q(tmp, complete = TRUE)[, set, drop = FALSE]
}

.DataPostProc <- function (DataObj, obj, loadings, scores, cl, bScores)
{
	idx <- order (obj, decreasing = TRUE)
	obj <- obj [idx]
	loadings <- loadings [,idx, drop = FALSE]

	if (bScores)
		scores <- scores [,idx, drop = FALSE]

	ret <- list()

   ##loadings
	{
		c <- ncol (loadings)
		r <- nrow (loadings)
		ret$loadings <- loadings

		ret$loadings <- .loadSgnU (ret$loadings)

		if (is.null (dimnames (DataObj$x)[[2]]))
			dimnames (ret$loadings) <- list (paste (rep ("V", r), 1:r, sep = ""), paste (rep ("Comp.", c), 1:c, sep = ""))
		else
			dimnames (ret$loadings) <- list (dimnames (DataObj$x)[[2]], paste (rep ("Comp.", c), 1:c, sep = ""))

		class (ret$loadings) <- "loadings"
	}

   ##sdev
	ret$sdev <- as.numeric (obj)
	names (ret$sdev) <- dimnames (ret$loadings)[[2]]

   ##center
	ret$center <- DataObj$center
   ##scale
	ret$scale <- DataObj$scale
   ##n.obs
	ret$n.obs <- nrow (DataObj$x)

   ##scores
	if (bScores)
	{
		ret$scores <- scores
		dimnames (ret$scores) <- list (1:nrow (scores), dimnames (ret$loadings)[[2]]) ;
	}
	else
		ret$scores <- NULL

	ret$call <- cl

	class (ret) <- c ("pcaPP", "princomp")
	return (ret)
}

