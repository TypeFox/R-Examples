
#cor.fk <- function (x, y = NULL, cor = TRUE)					#corissue
cor.fk <- function (x, y = NULL)
{
	cor = TRUE													#corissue
	if (is.null (y))
	{		##	x is expected to be matrix or data frame
		if (!is.matrix (x) && ! is.data.frame (x))
			stop ("x must be either numeric vector, matrix or data.frame.")
		{
			p <- ncol (x)
			dn <- colnames (x)
			ret <- diag (p)
			dimnames (ret) <- list (dn, dn)

			for (i in 1:p)
			{
				if (!cor)	##	calculating the diagonal elements. if cor == TRUE the diagonal elements are all equal to 1 and hence don't need to be calculated
#					ret[i, i] <- cor.fk (x[, i], x[, i], cor = cor)
					ret[i, i] <- cor.fk (x[, i], x[, i])	#corissue

				if (i == p)
					return (ret)

				ord <- order (x[, i])
				cur.x <- x[ord, i]
				for (j in (i+1):p)
					ret[i, j] <- ret[j, i] <- .cor.fk.2d (cur.x, x[ord,j], cor)
			}
		}
	}
	else
	{
		if (length (x) != length (y))
			stop ("x and y must have same length.")
		ord <- order (x)
		return (.cor.fk.2d (x[ord], y[ord], cor))
	}
}

.cor.fk.2d <- function (x, y, cor)
{
	if (length (x) != length (y))
		stop ("x and y must have same length.")

	ret <- .C ("kendallNlogN", PACKAGE = "pcaPP", NAOK = FALSE, DUP = TRUE		##	20130322 set DUP = TRUE - kendallNlogN implementation modifies x & y vectors!!
				, as.double (x)
				, as.double (y)
				, as.integer (c (length (x), cor))
				, ret = double (1)
				)
	return (ret$ret)
}
