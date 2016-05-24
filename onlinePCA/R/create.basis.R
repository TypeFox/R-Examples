create.basis <- function(x, p, sp = 1e-9, degree = 3, nderiv = 2)
{
 	d <- length(x)
	ord <- degree + 1
	knot <- quantile(x, seq_len(p-ord)/(p-degree))
	knot <- c(rep_len(min(x),ord), knot, rep_len(max(x),ord))

	# B-spline matrix
	B <- splines::splineDesign(knot, x, ord)
	
	# Weights for numerical integration (trapeze method)
	step <- diff(x)
	w <- (c(0,step) + c(step,0))/2
	
	# Gram matrix of B
	M <- crossprod(B, w * B)
	eigM <- eigen(M)
	vec <- eigM$vectors
	val <- eigM$values
	sqrtM <- vec %*% (t(vec) * sqrt(val))
	invsqrtM <- vec %*% (t(vec) / sqrt(val))
	
	# Penalty matrix
	if (sp>0)
	{
		DB <- splines::splineDesign(knot, x, ord, rep_len(nderiv, d))
		P <- crossprod(DB, w * DB)
	}
	
	# Matrix that maps functional data to their
	# (smoothed) coefficients in B-spline basis
	S <- if (sp>0) {tcrossprod(sqrtM %*% solve(M + sp * P), w * B) 
			} else {tcrossprod(invsqrtM, w * B)}

	list(B = B, S = S, invsqrtM = invsqrtM) 
}
