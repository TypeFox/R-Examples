.sPCAgrid.R <- function(x, k = 1, splitcircle = 10, maxiter = 10, method = mad, lambda = 0, cut.pc = TRUE, trace =  0, f.rho = .ident, glo.scatter = 0, ord.mod = 2, center, scale, pc.ini, k.ini, inc.scores = FALSE, HDred = c ("eigen", "qr", "svd", "svd.mean", FALSE))
{

# grid search for PCA, using Householder transformation for orthogonality of loadings
#
# x ... centered (!) data matrix
# k ... number of components to determine
# splitcircle ... number of directions to search
# maxiter ... maximum number of iterations
# method ... how to estimate the standard deviation
#        should be: "sd", "mad", "qn"
#
	stopifnot (is.numeric (trace))
	
	
	n <- nrow (x)
	p <- ncol (x)
	
	pHD <- 0	
	HDProj <- NULL

	f.HDred <- c ("eigen", "qr", "svd.mean", "svd", FALSE)
	HDred <- pmatch (HDred[1], f.HDred)
	if (is.na (HDred))
		HDred <- 0
	else
		HDred <- f.HDred[HDred]

    if(p > n)			## Dimension reduction for high dimensional datatsets
    {
	    if (HDred == "svd")
	    {
	        svdx <- svd(t(x))
	        x <- svdx$v %*% diag (svdx$d)
	        HDProj <- svdx$u
    	}
	    else if (HDred == "svd.mean")
	    {
		    x.m <- colMeans (x)
	        svdx <- svd(t(x) - x.m)
	        x <- svdx$v[, -n] %*% diag (svdx$d[-n])
	        HDProj <- svdx$u[, -n]
    	}
	    else if (HDred == "qr")
	    {
			qrtx <- qr (t (x))
			x <- t (qr.R (qrtx))
			HDProj <- qr.Q (qrtx)
		}
		else if (HDred == "eigen")
		{
			e <- eigen (cov (x))
			ev <- e$vectors[,order (e$values, decreasing = TRUE)[1:(n-1)]]
			x <- x %*% ev
			HDProj <- ev[]
		}
		
		if (!is.null (HDProj))
		{
			p <- ncol (x)
			n <- nrow (x)
#cat ("reduced dimensions -> n x p =", n, "x", p, "\n")
		}
   	}

   	stopifnot (k >= 1)

	len.lambda = length (lambda)

	if (len.lambda != 1 && len.lambda != k)
		warning ("length (lambda) should either be equal to 1 or k")
	lambda <- rep (lambda, len = k)

	if (!missing (pc.ini) && !is.null (pc.ini))
	{
		if (missing (k.ini))
			k.ini <- pc.ini$k
	}
	else
		k.ini <- 0

	if (k.ini)
	{
		stopifnot (nrow (pc.ini$load) == ncol (x))
		stopifnot (ncol (pc.ini$load) == ncol (x))

		if (missing (k.ini))
			k.ini <- pc.ini$k
		stopifnot (k.ini + k <= ncol (x))

		k <- k.ini + k
		
		stopifnot (length (pc.ini$sdev) == ncol (x))

		sdev <- pc.ini$sdev

		l <- pc.ini$load
	}
	else
	{
		k.ini <- 0
		sdev <- rep (NA, p)
		l <- diag (p)
	}

	if (missing (center))
		center = NULL
	if (missing (scale))
		scale = NULL

	x <- ScaleAdv (x, center, scale)$x

	stopifnot (k <= ncol (x))

# 	if (splitcircle %% 2 == 0)				##	forcing splitcircle to be odd (-> obj. function always increases during maxiter - apart from numerical issues)...
# 	#if (splitcircle %% 2)					##	forcing splitcircle to be even
# 		splitcircle = splitcircle + 1

	if (k.ini)
		y <- x %*% l[,(k.ini + 1):p]
	else
		y <- x
		
	
	if (glo.scatter == 0)
		scl.mean <- sqrt (mean (apply (x, 2, method)^2))
	else
		scl.mean <- 1

	for (nb in (k.ini + 1):k)
	{  # loop over number of comp
		if (glo.scatter == 1)
			scl.mean <- sqrt (mean (apply (y, 2, method)^2))

		p1 <- p - nb + 1 # dimension will be reduced for subsequent PCs			##	p1 = dimensionality of remaining subspace

		if (p1 == 1)	##	if the subspace is one dimensional the only thing left to do is to calculate sdev of the remaining (nx1) matrix Y
		{				##	the loadings are already fixed - there's no free parametes (d.o.f) left
			sdev [nb] = method (y)
			next (nb)
		}

			# ordering of variables by variance

		nord <- order(apply(y,2,method),decreasing=TRUE)				##	get the decreasing order of the scale estimates.

		yord <- y[,nord]												##	order the data regarding to their scale..

			#	initializing the lincomb
		afin <- c(1, rep(0,p1-1)) 										##	with (1, 0, 0, 0, ...)
		yopt = yord[,1]													##	y %*% afin

			# outer loop to run until convergence

		for (i in 0:maxiter)											##	the 0 round initializes the system (as before), with div =2^i = 1
		{
			sumabsdelta = 0
			for (j in 1:p1)
			{
				if (abs (afin [j]) == 1)	##  if abs (afin[j]) == 1 -> only the jth variable is considered so far (norm (afin) == 1)-> 
					next (j)				##		it wouldn't change anything adding the jth component again. we would only run into numerical issues.

				rf = c(sqrt (1 - afin[j]^2), afin[j])	## rotation factors = cos (alpha), sin (alpha)

				if (afin [j] != 0)			##	 setting afin[j] to zero. only necessary when afin [j] != 0
				{							##		-> changing afin[-j] accordingly
					yopt <- (yopt - yord[,j] * rf[2]) / rf[1]
					afin [-j] <- afin [-j] / rf[1]
					afin [j] <- 0
				}

				cury <- cbind(yopt,yord[,j])
#				if (glo.scatter >= 3)
#					scl.mean <- sqrt (mean (apply (cury, 2, method)^2))

				res <- .gridplane.shrk(n = splitcircle, div = 2^i, curL = rf[2], y = cury, method = method, afin = afin, nord = nord, curP = j, curK = nb, l = l, lambda = lambda[nb - k.ini], f.rho = f.rho, scl.mean = scl.mean, u = HDProj)

				sumabsdelta = sumabsdelta + abs (res$alphamax[2] - rf[2])	##	sums up the absolute change of the angles

				rf <- res$alphamax			## rotation factors <- c (cos (alpha), sin (alpha))

				yopt <- yopt * rf[1] + yord[,j] * rf[2]
				afin[-j] = afin[-j] * rf[1]
				afin[j] = rf[2]

#cat ("loadings:", afin, "\n")
			}

			afin=afin/sqrt(sum(afin^2))		##	normieren. -> why?
			objf <- res$objmax

			{		##	thinking about dropping this block. if splitcircle odd -> the obj can't drop (or?!?!?)! (because the best candidate in the last round would hence be a candidate in this round too)
				if (!i || objf>=objfold)
				{
					objfbest <- objf
					objfold <- objf
					sclbest <- 	res$sclmax
					afinbest <- afin/sqrt(sum(afin^2)) ##  why again?
				}
				else if (trace >= 2)
					cat ("objective function dropped:", objfold - objf, "\r\n")
			}
			if (sumabsdelta < 1e-16)
			#if (2^-i < 1e-16)
			{
				cat ("stopping execution after", i + 1, "loops\n")
				break ## breaks the maxiter loop, as there's no more progress..
			}
		} # maxiter

		sdev[nb] <- sclbest

# 		if (length (afinbest) == 2)
# 		{
# 			afin[nord] <- afinbest
# 			afinNV <- matrix (c (afin, afin[2], -afin[1]), nrow = 2)
# 		}
# 		else
		{
		    afinNV <- diag(p1)
	
		    if (ord.mod == 1)
		    {
		    	afin <- afinbest				##	xxxord1
		    	idx.ref <- 1					##	xxxord1
	    	}
		    else
		    {
		    	afin[nord] <- afinbest			##	xxxord2
		    	idx.ref <- nord[1]#which.max (abs (afin))##	xxxord2
	    	}

		    N <- afinNV[,idx.ref]-afin

		    N.norm <- .norm (N)

		    if (N.norm > 1e-6)		##	use zero.tol
		    {
				if (afin [idx.ref] < 0)	##	we want the first (one particular - which?) value to be positive
					afin <- -afin
				N <- N / N.norm

				afinNV <- afinNV - (2*N %*% t(N))
		    }

			#afinNV <- cbind (afin, Null (afin))		##	Null creates the orthogonal complement of afin	    

			if (ord.mod == 1)
			{
				afinNV[nord,] <- afinNV		##	xxxord1
														##	nord must be inversed HERE, not before creating the backtransformation.
														##	otherwise the result would depend on the input order of variables (due to numerical issues when creating the back transformation matrix)
														##		thus ord.mod == 1 is the "better" way? why is default 2 then? maybe because the "real" PCAgrid algo does it that way?
			}
			else
				afinNV <- afinNV[,nord]			##	xxxord2
		}

		## afinNV contains afin + its orthogonal complement
		l[,nb:p] <- l [,nb:p] %*% afinNV

		y <- y %*% afinNV[,-1]						##	projecting the data into the complement of afin.

		if (trace >= 3)
		{
			y1=x%*%l[,(nb+1):p]
			cat ("numerical instability for component", nb, ":", .norm (y-y1), "\r\n")
		}
	} # loop over number of comp

	ord <- order(sdev,decreasing=TRUE)
	sdev <- sdev[ord]

	l <- l [,ord, drop = FALSE]
	
	proj.loadings <- NULL
	if(!is.null (HDProj))			##	undo Projection for high dimensional datasets
	{
#		proj.loadings <- l
		l <- HDProj %*% l
	}

	if (cut.pc)
	{
		l <- l [,1:k, drop = FALSE]
		sdev <- sdev [1:k]
	}

	ret <- list (sdev = sdev, loadings = l, scores = NULL, k = k, call = match.call) #, u = HDProj, proj.loadings = proj.loadings, x = x)

	if (inc.scores)
		ret$scores <- x %*% l

	ret	
}

.gridplane.new <- function (n = 10, div = 1, curL = 0, y, method, ...)
{
	nangle <- seq (-pi/2, pi/2, len = n) / div	+	##	the angles to check
					asin (curL)						##	the angle - component of the current loading
	alpha <- cbind(cos(nangle),sin(nangle))

	obj <- apply(y %*% t(alpha), 2, method)											##	VVV like gridplane VVV
	idx.max = which.max (obj)
	list(objmax=obj[idx.max],alphamax=alpha[idx.max,], sclmax=obj[idx.max], idx.max = idx.max)
}

.ident <- function (i) 1
.inverse <- function (i) 1/(i+1)

.gridplane.shrk <- function (n = 10, div = 1, curL = 0, y, method, afin, nord, curP, curK, l, lambda, f.rho = .ident, scl.mean, u, ...)
{
	nangle <- seq (-pi/2, pi/2, len = n) / div	+	##	the angles to check
					asin (curL)						##	the angle - component of the current loading
	alpha <- cbind(cos(nangle),sin(nangle))

													##	afinN the new linear combinations, one column for each angle tested...
	afinN <- afin %*% t(alpha[,1])						##	applying the effect of the new param "curP" to all other params
	afinN[curP,] <- alpha[,2]							##	setting the new param "curP"

	afinN[nord,] <- afinN							##	reversing the sort order of the variables

	p <- nrow (l)

	
	#nOrigParams <- l [, curK:p] %*% afinN			##	transforming the new params back into the "real" system of Coordinates

	if (!missing (u) && !is.null (u))				##	undo the svd transformation for high dimensional datasets
		u <- u %*% l [, curK:p]
	else
		u <- l [, curK:p]
	nOrigParams <- u %*% afinN

# cat ("backtrans: \n")
# print (u)
# cat ("\n\n")
	
	shrk <- colSums (abs (nOrigParams))				##	the shrinkage factors for each angle

	scl <- apply(y %*% t(alpha), 2, method)			##	the scale for each angle
	obj <- scl^2 - lambda * shrk * scl.mean^2 * f.rho (curK)		##	objective function...
#	obj <- - shrk
#browser ()	
	idx.max = which.max (obj)
	nang <- nangle[idx.max]
# 	if (nang < 0.4)
# 		nang = - (nang - pi / 2)

#cat ("checking angles ", paste (nangle, " (", - lambda * shrk * scl.mean^2, ", ", obj, "),", collapse = " ", sep = ""), "\n", sep = "")



#cat ("selected angle", format (nang, digits = 22), "\n")
	list(objmax = obj[idx.max], alphamax = alpha [idx.max,], sclmax = scl[idx.max], idx.max = idx.max)
}

.norm <- function (x) sqrt (sum(x^2))
.sumsq <- function (x) sum (x^2)
