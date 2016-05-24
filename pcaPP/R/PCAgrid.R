
sPCAgrid <- function(x, k = 2, method = c ("mad", "sd", "qn"), lambda = 1#, norm.q = 1, norm.s = 1
					, maxiter = 10, splitcircle = 25, scores = TRUE, zero.tol = 1e-16, center = l1median, scale, trace = 0, store.call = TRUE, control, ...)
{

	norm.q <- .get_par (list (...), "norm.q", 1)
	norm.s <- .get_par (list (...), "norm.s", 1)
	glo.scatter <- .get_par (list (...), "glo.scatter", 0)

	check.orth <- FALSE
	dat <- list (x = x, substitute_x = substitute (x), k = k, method = method, maxiter = maxiter, splitcircle = splitcircle, check.orth = check.orth, scores = scores, lambda = lambda, zero.tol = zero.tol, center = center, store.call = store.call, trace = trace, ...)

	if (!missing (scale))	dat$scale <- scale

	if (!missing (control))
		dat <- .ParseControlStructureC (dat, control)

#	dat$check.orth <- FALSE		##	2do: remove this value!
	dat$HDred <- FALSE
	dat$cut.pc <- TRUE
    dat$glo.scatter <- glo.scatter
	dat$SpeedUp <- 0

	dat <- .sPCAgrid..DataPreProc (dat)
	dat$call <- match.call ()

	n <- nrow (dat$x)

	ret.C <- .C ("sPCAgrid", PACKAGE = "pcaPP", NAOK = TRUE
				, nParIn = as.integer (c(dim (dat$x), dat$k, dat$splitcircle, dat$maxiter, dat$method, dat$trace, dat$k.ini, dat$check.orth, dat$glo.scatter, dat$pHD, dat$SpeedUp))
				, nParOut = integer (1)
				, dParIn = as.double (c (dat$zero.tol, norm.q, norm.s))
				, as.double (dat$x)
				, l = as.double (dat$l)
				, sdev = as.double (dat$sdev)
				, obj = double (dat$k)
#				, max.maha = double (n)
				, as.double (dat$lambda)			## length = k - k.ini + 1
				, as.double (dat$HDProj)
				)

	return (.sPCAgrid..PostProc (dat, ret.C))
}

PCAgrid <- function (x, k = 2, method = c ("mad", "sd", "qn"), maxiter = 10, splitcircle = 25, scores = TRUE, zero.tol = 1e-16, center = l1median, scale, trace = 0, store.call = TRUE, control, ...)
{
	check.orth <- FALSE
	dat <- list (x = x, substitute_x = substitute (x), k = k, method = method, maxiter = maxiter, splitcircle = splitcircle, check.orth = check.orth, scores = scores, zero.tol = zero.tol, center = center, store.call = store.call, trace = trace, ...)

	if (!missing (scale))	dat$scale	<- scale

	if (!missing (control))
		dat <- .ParseControlStructureC (dat, control)

#	dat$check.orth <- FALSE		##	2do: remove this value!
	dat$HDred <- TRUE
	dat$cut.pc <- TRUE
	dat$glo.scatter <- 0
	dat$SpeedUp <- 0

	dat <- .PCAgrid..DataPreProc (dat)
	dat$call <- match.call ()

	n <- nrow (dat$x)

	ret.C <- .C ("PCAgrid", PACKAGE = "pcaPP", NAOK = TRUE
				, nParIn = as.integer (c(dim (dat$x), dat$k, dat$splitcircle, dat$maxiter, dat$method, dat$trace, dat$k.ini, dat$check.orth))
				, nParOut = integer (1)
				, dParIn = as.double (dat$zero.tol)
				, as.double (dat$x)
				, l = as.double (dat$l)
				, sdev = as.double (dat$sdev)
				, obj = double (dat$k)
#				, max.maha = double (n)
				)

	return (.PCAgrid..PostProc (dat, ret.C))
}

.sPCAgrid.ini <- function	(x, k = 2, method = c ("mad", "sd", "qn"), norm.q = 1, norm.s = 1, maxiter = 10, splitcircle = 25, scores = TRUE, zero.tol = 1e-16, center = l1median, scale, store.call = TRUE, trace = 0, cut.pc = TRUE, 
							pc.ini, k.ini, ord.all = FALSE, HDred = FALSE, lambda = 1, glo.scatter = 0, SpeedUp = 0, check.orth = FALSE, control, ...)
{
	dat <- list (x = x, substitute_x = substitute (x), k = k, method = method, maxiter = maxiter, splitcircle = splitcircle, scores = scores, lambda = lambda, zero.tol = zero.tol, center = center, store.call = store.call, trace = trace, cut.pc = cut.pc, glo.scatter = glo.scatter, ord.all = ord.all, HDred = HDred, SpeedUp = SpeedUp, check.orth = check.orth, ...)

	if (!missing (scale))	dat$scale <- scale
	if (!missing (pc.ini))	dat$pc.ini <- pc.ini
	if (!missing (k.ini))	dat$k.ini <- k.ini

	if (!missing (control))
		dat <- .ParseControlStructureC (dat, control)

	dat <- .sPCAgrid..DataPreProc (dat)
	dat$call <- match.call ()

	n <- nrow (dat$x)

	ret.C <- .C ("sPCAgrid", PACKAGE = "pcaPP", NAOK = TRUE
				, nParIn = as.integer (c(dim (dat$x), dat$k, dat$splitcircle, dat$maxiter, dat$method, dat$trace, dat$k.ini, dat$check.orth, dat$glo.scatter, dat$pHD, dat$SpeedUp))
				, nParOut = integer (1)
				, dParIn = as.double (c (dat$zero.tol, norm.q, norm.s))
				, as.double (dat$x)
				, l = as.double (dat$l)
				, sdev = as.double (dat$sdev)
				, obj = double (dat$k)
#				, max.maha = double (n)
				, as.double (dat$lambda)			## length = k - k.ini + 1
				, as.double (dat$HDProj)
				)

	return (.sPCAgrid..PostProc (dat, ret.C))
}

.ParseControlStructureC <- function (x, control)
{
	if (is.null (control) ||
		!length (control))
		return (x)

	nc <- names (control)
	if (is.null (nc) ||
		any (nc == ""))
		stop ("Each item of list \x22control\x22 must have a name.")

	for (i in 1:length (control))
		x[[nc[i]]] <- control[[i]]
	return (x)
}

.sPCAgrid..DataPreProc <- function (x)
{
	len.lambda = length (x$lambda)

	if (len.lambda != 1 && len.lambda != x$k)
		warning ("length (lambda) should either be equal to 1 or k")
	x$lambda <- rep (x$lambda, len = x$k)

	x <- .PCAgrid..DataPreProc (x)

	return (x)
}

.validScaleMethods <- function () c ("sd", "mad", "qn") 
.getScaleMethod <- function (x)
{
	.validMethods <- c ("sd", "mad", "qn")
	x <- x[1]
	if (is.character (x))
	{
		method <- match.arg (x, .validScaleMethods ())
		return (match (method, .validScaleMethods()) - 1)
	}

	if (x >= 3 
		&& x != 5                                                               ## hack for the sPCAgrid-paper -> remove this line again...
		) 
		stop ("the method is supposed to be a value < 3")

	return (x)
}

.getScaleName <- function (x)
{
	if (is.function (x))
		return (.GetFunctionName (x))

	return (.validScaleMethods () [.getScaleMethod (x) + 1])
}


.getScaleFunction <- function (x)	##	was .EvalScaleFunction before
{
	if (is.function (x))
		return (x) 

	f.idx <- .getScaleMethod (x) + 1

	F <- c (sd, mad, qn)
	return (F [[f.idx]])
#	return  (eval (parse (text = method)))
}

.Check_DimRed <- function (x)
{
		##	dimension reduction
	x$pHD <- 0			##	stores the original p, if we have to reduce the data matrix' dimensionality (n < p)

	x$x.orig <- x$x

    if(x$p > x$n && x$HDred)			## Dimension reduction for high dimensional datatsets
    {
		svdx <- svd(t(x$x))
		x$svdx <- svdx

		x$x <- svdx$v %*% diag (svdx$d)
		x$HDProj <- svdx$u

		x$pHD <- x$p
		x$p <- ncol (x$x)
		x$n <- nrow (x$x)
		if (x$trace >= 2)
			cat ("reduced dimensions -> n x p =", x$n, "x", x$p, "\n")
   	}

	return (x)
}


.scale <- function(x)
{
	x$scl <- ScaleAdv (x$x, x$center, x$scale)
	x$x <- x$scl$x

	if (x$pHD) # center and scale must have original data dimension:
	{
		x$scl$center <- as.vector(x$svdx$u%*%x$scl$center)
		x$scl$scale <- ScaleAdv(x$x%*%t(x$svdx$u), center = NULL, scale = x$scale)$scale
	}

	return (x)
}


.Check_pc.ini <-function (x)
{
	if (!is.null (x$pc.ini))
	{
		if (is.null (x$k.ini))
			x$k.ini <- x$pc.ini$k
	}
	else
		x$k.ini <- 0

	if (x$k.ini)
	{
		stopifnot (all (dim (x$pc.ini$load) == ncol (x$x)))
		stopifnot (x$k.ini + x$k <= ncol (x$x))
		if (is.null (x$k.ini))
			x$k.ini <- x$pc.ini$k

		x$k <- x$k.ini + x$k	##	2do: don't change k anymore to k.ini + k
		stopifnot (length (x$pc.ini$sdev) == ncol (x$x))
		x$sdev <- x$pc.ini$sdev
		x$l <- x$pc.ini$load
	}
	else
	{
		x$k.ini <- 0
		x$sdev <- rep (NA, x$p)
		x$l <- diag (x$p)
	}

	return (x)
}

.PCAgrid..DataPreProc <- function (x)
{

	x$args <- list (method = x$method)		## 2do: which further args go here?
	##x[names (x) != "x"]	##	storing the function arguments

	x$x <- X <- .Conv2Matrix (x$x, x$substitute_x)

	if (is.function (x$method) && !is.null (x$method.name))
		attributes (x$method)$NAME <- x$method.name

	##	checking/initializing parameters

	x$n <- nrow (x$x)
	x$p <- ncol (x$x)

	stopifnot (is.numeric (x$trace))
	stopifnot (x$k >= 1)
	stopifnot (length (x$k) == 1)

	x$method <- .getScaleMethod (x$method)

	x <- .Check_DimRed (x)	#
	x <- .Check_pc.ini (x)

	stopifnot (x$k <= ncol (x$x))

	x <- .scale (x)

	return (x)
}

.sPCAgrid..GetLambda.ini <- function (x)
{
	if (!x$k.ini)
		return (NULL)
	if (is.null (x$pc.ini$lambda))
		return (rep (0, x$k.ini))
	return (x$pc.ini$lambda)
}

.sPCAgrid..PostProc <- function (x, ret.C)
{
	ret <- .PCAgrid..PostProc (x, ret.C)

	ret$lambda <- c (.sPCAgrid..GetLambda.ini (x), x$lambda)

	ret
}

.cut.pc <- function (ret, k)
{
	ret$loadings <- ret$loadings [, 1:k, drop = FALSE]
	ret$sdev <- ret$sdev [1:k]
	if (!is.null (ret$scores))
		ret$scores <- ret$scores[, 1:k, drop = FALSE]
	if (!is.null (ret$lambda))
		ret$lambda <- ret$lambda [1:k]

	return (ret)
}

.orderPCs <- function (ret, k, k.ini, ord.all = FALSE)
{
	if (is.null (ord.all) || ord.all)
		idx.ord <- 1:k
	else
		idx.ord <- (k.ini + 1):k

	if (length (idx.ord) == 1)	##	nothing to sort.
		return (ret)

	ord <- order(ret$sdev[idx.ord], decreasing = TRUE)
	ord <- idx.ord[ord]

	ret$pc.order <- 1:k
	ret$pc.order[idx.ord] <- ord

	ret$sdev[idx.ord]  <- ret$sdev[ord]
	ret$loadings[, idx.ord]  <- ret$loadings [, ord]
	ret$obj[idx.ord]  <- ret$obj [ord]

	if (!is.null (ret$lambda))
		ret$lambda[idx.ord]  <- ret$lambda [ord]


	return (ret)
}

.PCAgrid..PostProc <- function (x, ret.C)
{

	ret <- list (sdev = ret.C$sdev, 
				loadings = matrix (ret.C$l, ncol = x$p, nrow = x$p), k = x$k, 
				obj = ret.C$obj, 
				n.obs = nrow (x$x.orig), args = x$args, #call = x$call,
				scale = x$scl$scale, center = x$scl$center
#				, max.maha = ret.C$max.maha
				)

	if (x$store.call)
		ret$call = x$call

	if(x$pHD)			##	undo SVD for high dimensional datasets
		ret$loadings <- x$HDProj %*% ret$loadings

	if (x$cut.pc)
		ret <- .cut.pc (ret, x$k)

	ret$loadings <- .loadSgnU (ret$loadings)

	ret <- .orderPCs (ret, x$k, x$k.ini, x$ord.all)

	ndn <- list (NULL, paste ("Comp", 1:ncol(ret$loadings), sep = "."))	##	new dimnames

	if (is.null (dimnames (x$x)[[2]]))
		ndn[[1]] <- paste ("X", 1:nrow(ret$loadings))
	else
		ndn[[1]] <- dimnames (x$x.orig)[[2]]

	dimnames (ret$loadings) <- ndn

	if (x$scores)
	{
		ret$scores <- t (t (x$x.orig) - x$scl$center) %*% ret$loadings
		dimnames (ret$scores) <- list (dimnames (x$x)[[1]], ndn[[2]])
	}

	class (ret$loadings) <- "loadings"
	class (ret) <- "princomp"

	return (ret)
}

.get_par <- function (l, idx, default)
{
	ret <- l[[idx]]
	if (is.null (ret))
		return (default)
	return (ret)
}

#.HDred <- function ()
#{
#	if (HDred == "svd")
#	{
#		svdx <- svd(t(x))
#		x <- svdx$v %*% diag (svdx$d)
#		HDProj <- svdx$u
#	}
#	else if (HDred == "svd.mean")
#	{
#	x.m <- colMeans (x)
#	svdx <- svd(t(x) - x.m)
#	x <- svdx$v[, -n] %*% diag (svdx$d[-n])
#	HDProj <- svdx$u[, -n]
#	}
#	else if (HDred == "qr")
#	{
#	qrtx <- qr (t (x))
#	x <- t (qr.R (qrtx))
#	HDProj <- qr.Q (qrtx)
#	}
#	else if (HDred == "eigen")
#	{
#	e <- eigen (cov (x))
#	ev <- e$vectors[,order (e$values, decreasing = TRUE)[1:(n-1)]]
#	x <- x %*% ev
#	HDProj <- ev
#	}
#
#	if (!is.null (HDProj))
#	{
#	pHD <- p
#	p <- ncol (x)
#	n <- nrow (x)
#	if (trace >= 2)
#		cat ("reduced dimensions -> n x p =", n, "x", p, "\n")
#	}
#}
