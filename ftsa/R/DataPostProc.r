DataPostProc <- function (DataObj, obj, loadings, scores, cl, bScores)
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

