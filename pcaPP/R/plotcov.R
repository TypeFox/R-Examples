plotcov <- function (cov1, cov2, method1, labels1, method2, labels2, ndigits = 4, ...)
{

	if (class (cov1) == "matrix")
		cm1 = cov1
	else if (is.null (cov1$cov))
		stop ("No appropriate covariance structure specified")
	else
	{
		cm1 = cov1$cov
		if (!is.null (cov1$method))
			method1 = cov1$method
	}

	if (missing (method1)) #|| is.null(method1))
		method1 = "Method 1"

	if (ncol (cm1) != nrow (cm1))
		stop ("Supplied covariance structure has to be quadratic!")

	if (missing (labels1))
		if (is.null (dimnames (cm1)))
			labels1 = paste ("V", 1:ncol (cm1), sep = "")
		else
			labels1 = dimnames (cm1)[[1]]

	if (missing (cov2))
	{	# only plotting cov1
		.plotsingle (cm1, labels1, method1, ndigits, ...)
	}
	else
	{
		if (class (cov2) == "matrix")
			cm2 = cov2
		else if (is.null (cov2$cov))
			stop ("No appropriate covariance structure specified")
		else
		{
			cm2 = cov2$cov
			if (!is.null (cov2$method))
				method2 = cov2$method
		}
		if (missing (method2))# || is.null(method2))
			method2 = "Method 2"

		if (ncol (cm2) != nrow (cm2))
			stop ("Supplied covariance structure has to be quadratic!")

		if (missing (labels2))
			if (is.null (dimnames (cm2)))
				labels2 = paste ("V", 1:ncol (cm2), sep = "")
			else
				labels2 = dimnames (cm2)[[1]]

		.plotcomp (cm1, cm2, labels1, labels2, method1, method2, ndigits, col = "#0096FF", ...)
	}
	invisible()
}

.plotcomp <- function (cor1, cor2, labels1, labels2, method1, method2, ndigits = 4, ...)	## internal function
{
	plot.new()
# 	plot (0,0, pch = "", axes = F, xlab = "", ylab = "")
	oldmar = par ("mar")
	par (mar = rep(1, 4))
	p = ncol (cor1)
	lim = c(-1, p + 1)
	plot.window(xlim=lim, ylim= lim, xaxs="i", yaxs="i")

	for (i in 1:p)
	{
		text (i - 0.5, p + 0.5, labels1[i], srt=90)
		text (-0.5, p - i + 0.5, labels2[i])
	}

	for (i in 2:p)
	{
		for (j in 1:(i-1))
		{
			.doEllipses (cor1[c(i,j), c(i,j)], pos = c(i - 1.5,p - j - 0.5), lwd = 2)
			.doEllipses (cor2[c(i,j), c(i,j)], pos = c(i - 1.5, p - j - 0.5), lwd = 2, ...)

			text (j - 0.5,p - i + 0.5, round (cor1[i,j], ndigits), adj = c(0.5,-0.1))
			text (j - 0.5,p - i + 0.5, round (cor2[i,j], ndigits), adj = c(0.5,1.1), ...)
		}
	}

	lines (c(0.5, p-0.5), c(p - 0.5, 0.5), lwd = 3)

#	lines (mean (lim) + 

	lines (lim[2] - 2 + c(-0.5, 0.3), c(-0.3, -0.3))
	lines (lim[2] - 2 + c(-0.5, 0.3), c(-0.7, -0.7), ...)
	text (lim[2] - 2 - 0.7, -0.275, method1, pos = 2)
	text (lim[2] - 2 - 0.7, -0.675, method2, pos = 2)

	par (mar = oldmar)
}


.plotsingle <- function (cm, labels, method, ndigits, ...)
{
	plot.new()

	oldmar = par ("mar")

	par (mar = rep(1, 4))
	p = ncol (cm)
	lim = c(-1, p + 1)

	plot.window(xlim = lim, ylim = lim, xaxs = "i", yaxs = "i")

	for (i in 1:p)
	{
		text (i - 0.5, p + 0.5, labels[i], srt=90)
		text (-0.5, p - i + 0.5, labels[i])
	}

	for (i in 2:p)
	{
		for (j in 1:(i-1))
		{
			.doEllipses (cm[c(i,j), c(i,j)], pos = c(i - 1.5,p - j - 0.5), lwd = 2, ...)
			text (j - 0.5 ,p - i + 0.5, round (cm[i,j], ndigits), ...)
		}
	}

	lines (c(0.5, p-0.5), c(p - 0.5, 0.5), lwd = 3)

	lines (lim[2] - 2 + c(-0.5, 0.3), c(-0.5, -0.5), col = col, ...)
	text (lim[2] - 2 - 0.7, -0.475, method, pos = 2)

	par (mar = oldmar)
}


.doEllipses <- function (acov, pos, ...)		## internal function
{
	acov = cov2cor (acov)
	cov.svd <- svd(acov, nv = 0)
	r <- cov.svd[["u"]] %*% diag(sqrt(cov.svd[["d"]]))

	m <- 100

	alphamd <- c(1/3)

#	par(pty="s")

    e1md <- cos(c(0:m)/m * 2 * pi) * alphamd
    e2md <- sin(c(0:m)/m * 2 * pi) * alphamd
    emd <- cbind(e1md, e2md)
    ttmd <- t(r %*% t(emd)) + rep(1, m + 1)
    lines(ttmd[, 1] + pos[1], ttmd[, 2]+ pos[2], ...)
}
