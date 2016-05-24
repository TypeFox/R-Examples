
#-----------------------------------------------------------------------------------------
#  distance, with expanded functionality.
#
#  "biom" method is deliberately prototyped with only the arguments
#  that it actually touches, plus those (i.e., "method") preceding.
#-----------------------------------------------------------------------------------------

distx <- function (x, ...) UseMethod ("distx")

distx.biom <- function (x, method="euclidean", groups=NULL, ..., bycol=TRUE) {
	distx(
		as.matrix (x, expand=TRUE), 
		method, 
		if (bycol) subColumn (groups, x) else subRow (groups, x),
		...,
		bycol=bycol)
	}

distx.matrix <- function(
	x, 
	method=c("euclidean", "bray-curtis", "jaccard", "mahalanobis", "sorensen", "difference", "maximum", "manhattan", "canberra", "binary", "minkowski"),
	groups=NULL, 
	p=NULL,
	..., 
	bycol=TRUE) {

	method <- match.arg (method)
	if (bycol) x <- t(x)
	dist.fun <- if (method %in% c ("bray-curtis", "jaccard", "mahalanobis", "sorensen", "difference")) {
		ecodist::distance
	} else stats::dist

	if (!is.null (p)) {
		dist2p <- apply (x, 1, 
				function (r, p, m) dist.fun (rbind (r, p), m),
				p, method)
		return(
			if (is.null (groups)) {
				dist2p
			} else tapply (dist2p, groups, mean))
		}

	if (is.null (groups)) return (dist.fun (x, method))

	groups <- as.factor (groups)
	D <- dist.fun (x, method)
	from <- unlist (sapply (2:nrow(x), seq, to=nrow(x)))
	to <- unlist (mapply (rep, 1 : (nrow(x)-1), (nrow(x)-1) : 1))
	zz <- tapply (D, list (groups [from], groups [to]), mean)

####  make symmetric by replacing NA's across the diagonal

	dd <- diag (zz)
	zz [is.na (zz)] <- 0
	zz <- zz + t(zz)
	diag (zz) <- dd
	zz
	}
