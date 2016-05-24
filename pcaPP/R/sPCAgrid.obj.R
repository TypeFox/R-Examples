####################################
##	Tradeoff Area - Optimization  ##
####################################

.TPO <- function (f.x = .l0sparse, f.y = .sumVar, obj.pc.0, obj.pc.1, fa.Name = "TPO", ...)
{
	if (missing (obj.pc.0))
		return (list (x = f.x (...), y = f.y (...)))
	if (missing (obj.pc.1))
		return (list (x = f.x (v.0 = obj.pc.0[[1]], ...), y = f.y (v.0 = obj.pc.0[[1]], ...)))

	x <- f.x (v.0 = obj.pc.0$x, v.1 = obj.pc.1$x, ...)
	y <- f.y (v.0 = obj.pc.0$y, v.1 = obj.pc.1$y, ...)

	return (-x * y)
}

.l1sparse <- function (x, pc, k, fa.Name = expression (paste (L[1], " Sparseness")), ...)
{
	p <- nrow (pc$load)
	n.k <- length (k)

	n.k * sqrt (p) - sum (abs (pc$loadings[,k]))
}

.l0sparse <- function (x, pc, k, zero.tol = 1e-16, fa.Name = expression (paste (L[0], " Sparseness")), ...)
{
	sum (abs (pc$loadings[,k]) <= zero.tol)
}

.pl1sparse <- function (x, pc, k, fa.Name = expression (paste (L[1], " Sparseness (%)")), ...)
{
	p <- nrow (pc$load)
	n.k <- length (k)

	.l1sparse (x, pc, k, ...) / ((sqrt (p) - 1) * n.k) * 100
}

.pl0sparse <- function (x, pc, k, zero.tol = 1e-16, fa.Name = expression (paste (L[0], " Sparseness (%)")), ...)
{
	p <- nrow (pc$load)
	n.k <- length (k)

	.l0sparse (x, pc, k, ...) / ((p - 1) * n.k) * 100
}


.lambda <- function (x, pc, k, fa.Name = expression (lambda), ...)
{
	if (length (k) == 1)
		return (pc$lambda[k])

	if (any (diff (pc$lambda) != 0) && length (k) > 1)
		stop ("This option can only be used when all considered PCs were calculated with the same lambda.")

	return (pc$lambda[1])
}

.vn <- function (n, k) 
{
	if (length (k) == 1)
		return (n)
	return (paste ("Cumulated", n))
}

.sumVar <- function (x, pc, k, NAME = FALSE, ...)
{
	if (NAME)
		return (.vn ("Explained Variance", k))
			
	return (sum (pc$sdev[k]^2))
}

.sumVarDiff <- function (x, pc, k, v.1 = 0, NAME = FALSE, ...)
{
	if (NAME)
		return (.vn ("diff Explained Variance", k))

	return (sum (pc$sdev[k]^2) - v.1)
}

.sumVarP <- function (x, pc, k, v.0, NAME = FALSE, ...)
{
	if (NAME)
		return (.vn ("Explained Variance (%)", k))

	if (missing (v.0))			##	return the total variance of x 
		return (sum (apply (x, 2, .getScaleFunction (pc$args$method))^2))

	return (sum (pc$sdev[k]^2) / v.0 * 100)
}

.sumVarDiffP <- function (x, pc, k, v.0, v.1 = 0, NAME = FALSE, ...)
{
	if (NAME)
		return (.vn ("diff Explained Variance (%)", k))

	if (missing (v.0))			##	return the total variance of x 
		return (sum (apply (x, 2, .getScaleFunction (pc$args$method))^2))

	return (sum (pc$sdev[k]^2) / v.0 * 100  - v.1)
}

.obj <- function (pc, fa.Name = "Objective Function", ...)
{
	pc$ev
}

#####################
##	BIC criterion  ##
#####################

.BIC.RSS <- function (obj.pc.0, f.BIC = .nBIC.Guo, f.RSS = .Calc.RSS.a, fa.Name = "BIC", ...)
#.BIC.RSS <- function (obj.pc.0, f.BIC = .nBIC.CC, f.RSS = .Calc.RSS.a, NAME = FALSE, ...)
#BIC.RSS <- function (obj.pc.0, f.BIC = .nBIC.CC, f.RSS = .Calc.RSS.e, get.name = FALSE, ...)
{
#	if (NAME)
#		return (paste (.GetFunctionName (f.BIC, ...), "(", .GetFunctionName (f.RSS, ...),")"))

	RSS <- f.RSS (...)

	if (missing (obj.pc.0))
		return (RSS)

	df_k <- .Calc.df_k (...)

	f.BIC (RSS0 = RSS, RSS = obj.pc.0, df_k = df_k, ...)
}

.Calc.df_k <- function (pc, k, zero.tol = 1e-15, ...) {	sum (abs (pc$loadings[,k]) > zero.tol)}

.Calc.RSS.e <- function (pc, k, fa.Name = "RSS: exact", ...)
{
	p <- nrow (pc$load)
	if (length (pc$sdev) < p)
		stop ("all eigenvectors needed for exact computation of RSS")

	sum (pc$sdev[-(1:max(k))]^2)	
}

.Calc.RSS.a <- function (x, pc, k, fa.Name = "RSS: approx", ...)
{
	if (ncol (pc$loadings) < length (k))
		stop ("at least k eigenvectors needed for approximate calculation of RSS")
	l <- pc$loadings[,1:max(k)]
	f.method <- .getScaleFunction (pc$args$method)
	sum (apply (x - x %*% l %*% t(l), 2, f.method)^2)
}

.nBIC.CC <- function (x, RSS0, RSS, df_k, k, fa.Name = "BIC", trace = 0, ...)
{
##	RSS0 as the restricted model's RSS
##	RSS as the unrestricted model's RSS
##	df_k as the number of non-zero loadings of the restricted model

	n <- nrow (x)
	BIC <- RSS0 / RSS + df_k * log (n) / n
	if (trace >= 2)
		cat ("RSS0:", RSS0, "; RSS:", RSS, "; df_k:", df_k, "; n:", n, "; BIC:", BIC, "\n")
	BIC	
}

.nBIC.Guo <- function (x, RSS0, RSS, df_k, k, fa.Name = "BIC", trace = 0, ...)
{
##	RSS0 as the restricted model's RSS
##	RSS as the unrestricted model's RSS
##	df_k as the number of non-zero loadings of the restricted model

	n <- nrow (x)
	BIC <- RSS0 / RSS + df_k * log (n)# / n
	if (trace >= 2)
		cat ("RSS0:", RSS0, "; RSS:", RSS, "; df_k:", df_k, "; n:", n, "; BIC:", BIC, "\n")
	BIC	
}
