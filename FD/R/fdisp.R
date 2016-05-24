`fdisp` <-
    function(d, a, tol = 1e-07)
{

if(!inherits(d, "dist"))
 stop("'d' must be a 'dist' object.")
n <- attr(d, "Size")

if (is.null(attr(d, "Labels") ) ) stop("'d' must have labels.","\n") else sn.d <- attr(d, "Labels")
if (missing(a) ){
  ab.names <- list("Community1", sn.d)
  a <- matrix(1, 1, n, dimnames = ab.names)
}
com <- nrow(a)
if (!is.matrix(a)) stop("'a' must be a matrix.")
if (ncol(a) != n )
 stop("Number of columns in 'a' must be equal to the number of objects in 'd'.")
if (is.null(colnames(a) ) ) stop("'a' must have column names", "\n") else sn.a <- colnames(a)
# check if species labels in 'd' and 'a' match
if (any(sn.d != sn.a) ) stop("Species labels in 'd' and 'a' need to be identical and ordered alphabetically (or simply in the same order).","\n")
# replace NA in 'a' by 0
a[which(is.na(a) )] <- 0
# check if one community has total abundance of zero (no species)
abun.sum <- apply(a, 1, sum)
if (any(abun.sum == 0) ) stop("At least one community has zero-sum abundances (no species).","\n")
# check if one species has total abundance of zero
abun.sum2 <- apply(a, 2, sum)
if (any(abun.sum2 == 0) ) stop("At least one species does not occur in any community (zero total abundance across all communities).","\n")
if (any(is.na(d) ) ) stop("NA's in the distance matrix.","\n")
A <- matrix(0, ncol = n, nrow = n)

A[row(A) > col(A)] <- -0.5 * d^2

A <- A + t(A)

# gower's double-centering
G <- bicenter.wt(A)
e <- eigen(G, symmetric = TRUE)

vectors <- e$vectors

eig <- e$values

# check if 'd' is Euclidean or not
w0 <- eig[n] / eig[1]
if (w0 > -tol) r <- sum(eig > (eig[1] * tol)) else r <- length(eig)
vectors <- vectors[, 1:r, drop = FALSE] %*% diag(sqrt(abs(eig <- eig[1:r])), r)

dimnames(vectors) <- list(colnames(a), NULL)
pos <- eig > 0

avg.dist.cent <- rep(NA, nrow(a)) ; names(avg.dist.cent) <- row.names(a)
for (i in 1:com){
  pres <- which(a[i ,] > 0)
  nb.sp <- nrow((unique(vec <- vectors[pres, , drop = F] ) ) )
  if (nb.sp >= 2){
    w <- a[i, pres]
    centroid <- apply(vec, 2, weighted.mean, w = w)
    dist.pos <- sweep(vec[, pos , drop = F], 2, centroid[pos])
    dist.pos <- rowSums(dist.pos^2)
    if (any(!pos) ){
	dist.neg <- sweep(vec[, !pos , drop = F], 2, centroid[!pos])
	dist.neg <- rowSums(dist.neg^2)
    }
    else dist.neg <- 0
    zij <- sqrt(abs(dist.pos - dist.neg) )
    avg.dist.cent[i] <- weighted.mean(zij, w)
  }
  else avg.dist.cent[i] <- 0
}
return(list(FDis = avg.dist.cent, eig = eig, vectors = vectors))
}
